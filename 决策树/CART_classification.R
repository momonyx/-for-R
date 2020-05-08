#定义Gini系数计算方式
giniy=function(y){
  return(1-sum((table(y)/length(y))^2))
}
ginixy=function(x,y,a){
  index=(x==a)
  return(giniy(y[index])*sum(index)/length(y)+giniy(y[!index])*sum(!index)/length(y))
}

#定义CART分类树节点形式
treenode=function(x,No=NA,father=NA,left=NULL,right=NULL,variable=NA,divchoice=NA,value=NA){
  return(list(No=No,father=father,left=left,right=right,variable=variable,divchoice=divchoice,value=value,index=x))
}

#CART分类树生成算法
#寻找使得gini系数最优特征和最优切分点
#叶节点条件：
#1.节点中的样本个数小于等于阈值(neps)
#2.gini(y)小于阈值(eps)->样本基本属于同一类
#3.没有更多可供选择的特征
cart.c=function(x,y,eps=1e-6,neps=1,cut=FALSE,alphaeps=NULL){
  n=nrow(x)
  p=ncol(x)
  xlevels=lapply(x,levels)
  
  #求解xi的使得Gini系数最小的切分点
  min.gini=function(x,y){
    xlevels=levels(x)
    gini=rep(NA,length(xlevels)) 
    for (i in 1:length(xlevels)){
      gini[i]=ginixy(x,y,xlevels[i])
    }
    gini[is.na(gini)]=Inf
    i=which.min(gini)
    return(c(gini[i],i))
  }
  
  currentnode=1
  tree=list()
  tree[[currentnode]]=treenode(1:n,1)
  
  while (TRUE){
    #更新当前节点所对应的类标记
    tree[[currentnode]]$value=names(which.max(table(y[tree[[currentnode]]$index])))
    
    #当节点中的样本个数小于阈值(neps)或gini(y)小于阈值(eps)时，将当前节点设置为叶节点
    if (giniy(y[tree[[currentnode]]$index])<eps || length(tree[[currentnode]]$index)<=neps){
      tree[[currentnode]]$left=NA
      tree[[currentnode]]$right=NA
      
    #其他情况需要计算gini(y|x)
    }else{
      
      #查找使得gini系数最小的最优特征和最优切分点
      gini=as.data.frame(lapply(x[tree[[currentnode]]$index,],min.gini,y[tree[[currentnode]]$index]))
      variable=which.min(gini[1,])
      
      #当没有更多可供选择的特征时，giniy(y[!index])=NaN->Inf，
      #所以若min(gini(y|x))=Inf则代表没有更多可供选择的特征，此时设置当前节点为叶节点
      if (gini[1,variable]==Inf){
        tree[[currentnode]]$left=NA
        tree[[currentnode]]$right=NA
        
      #其他情况当前节点均非叶节点
      }else{
        
        #更新最优特征和最优切分点
        tree[[currentnode]]$variable=variable
        tree[[currentnode]]$divchoice=xlevels[[variable]][gini[2,variable]]
        
        #更新当前节点的左右子节点，新建子节点
        #左子节点为满足分支条件的数据,右子节点为不满足分支条件的数据
        tree[[currentnode]]$left=length(tree)+1
        tree[[length(tree)+1]]=treenode(tree[[currentnode]]$index[x[tree[[currentnode]]$index,tree[[currentnode]]$variable]==tree[[currentnode]]$divchoice],
                                        No=length(tree)+1,father=currentnode)
        tree[[currentnode]]$right=length(tree)+1
        tree[[length(tree)+1]]=treenode(tree[[currentnode]]$index[x[tree[[currentnode]]$index,tree[[currentnode]]$variable]!=tree[[currentnode]]$divchoice],
                                        No=length(tree)+1,father=currentnode)
      }
    }
    
    #遍历树查询树找到下一更新节点
    for (i in 1:length(tree)){
      if (is.na(tree[[i]]$value)){newnode=tree[[i]]$No;break}
    }
    
    #当所有节点都进行过更新之后，跳出循环
    if (newnode!=currentnode) currentnode=newnode else break
  }
  if(!cut) return(tree) else return(cart.c.cut(x,y,tree,alphaeps))
}

#CART分类树剪枝算法
cart.c.cut=function(x,y,tree,alphaeps){
  alpha=rep(NA,length(tree))
  #对每一个子树计算满足L(father)=L(tree)的alpha
  for (i in 1:length(tree)){
    if (!is.na(tree[[i]]$left)){
      alpha[i]=giniy(y[tree[[i]]$index])-giniy(y[tree[[tree[[i]]$left]]$index])-giniy(y[tree[[tree[[i]]$right]]$index])
    }
  }
  #对于alpha大于阈值的子树，将其父节点设置为叶节点
  dellist=NULL
  for (i in which(alpha>alphaeps)){
    dellist=c(dellist,tree[[i]]$left,tree[[i]]$right)
    tree[[i]]$left=NA
    tree[[i]]$right=NA
  }
  #删除无用的子树
  if (!is.null(dellist)){
    for (i in dellist){
      tree[[i]]=0
    }
  }
  return(tree)
}


#CART分类树预测
cart.c.predict=function(xpre,tree){
  
  #对单条记录进行预测
  cart.c.predictONE=function(xpre,tree){
    #从根节点开始查找CART分类树
    currentnode=1
    #当其存在子节点时，若满足分支条件，则将当前节点更新左子节点，否则更新为右子节点
    while(!is.na(tree[[currentnode]]$left)){
      if (xpre[tree[[currentnode]]$variable]==tree[[currentnode]]$divchoice){
        currentnode=tree[[currentnode]]$left
      }else{currentnode=tree[[currentnode]]$right}
    }
    return(tree[[currentnode]]$value)
  }
  
  #对xpre用apply
  if (nrow(xpre)==1){return(cart.c.predictONE(xpre,tree))}else{
    return(apply(xpre,1,cart.c.predictONE,tree))
  }
}


#test part
x=matrix(c(c(rep("青年",5),rep("中年",5),rep("老年",5)),
           c("否","否","是","是","否","否","否","是","否","否","否","否","是","是","否"),
           c("否","否","否","是","否","否","否","是","是","是","是","是","否","否","否"),
           c("一般","好","好","一般","一般","一般","好","好","非常好","非常好","非常好","好","好","非常好","一般")),
         nrow=15)
x=as.data.frame(x)
y=c("否","否","是","是","否","否","否","是","是","是","是","是","是","是","否")
y=as.factor(y)

tree=cart.c(x,y)
treecut=cart.c.cut(x,y,tree,0.4)
cart.c(x,y,cut = TRUE,alphaeps = 0.4)

xpre=x[1,]
cart.c.predict(xpre,tree)
cart.c.predict(x,tree)
cart.c.predict(xpre,treecut)
cart.c.predict(x,treecut)
