#定义CART回归树节点形式
treenode=function(x,No=NA,father=NA,left=NULL,right=NULL,variable=NA,divchoice=NA,value=NA){
  return(list(No=No,father=father,left=left,right=right,variable=variable,divchoice=divchoice,value=value,index=x))
}

#最小二乘回归树生成算法
#寻找使得损失函数最优特征和最优切分点
#叶节点条件：
#1.节点中的样本个数小于等于阈值(neps)
#2.或(y-mean(y))^2小于阈值(eps)->样本基本属于同一类
#3.没有更多可供选择的特征
cart.r=function(x,y,eps=1e-6,neps=1,cut=FALSE,alphaeps=NULL){
  n=nrow(x)
  p=ncol(x)
  
  #回归树的损失函数
  lossfunc=function(x,y){
    minloss=Inf
    minA=0
    index=order(x)
    n=length(index)
    for (i in 1:(n-1)){
      if (all(!x[index[(i+1):n]]==x[index[i]])){
        c1=mean(y[index[1:i]])
        c2=mean(y[index[(i+1):n]])
        loss=sum((y[index[1:i]]-c1)^2)+sum((y[index[(i+1):n]]-c2)^2)
        if (loss<minloss) {minloss=loss;minA=x[index[i]]}
      }
    }
    return(c(minloss,minA))
  }
  
  currentnode=1
  tree=list()
  tree[[currentnode]]=treenode(1:n,1)
  
  while (TRUE){
    #更新当前节点所对应的类标记
    tree[[currentnode]]$value=mean(y[tree[[currentnode]]$index])
    
    #当节点中的样本个数小于阈值(neps)或(y-mean(y))^2小于阈值(eps)时，将当前节点设置为叶节点
    if ((y[tree[[currentnode]]$index]-tree[[currentnode]]$value)^2<eps || length(tree[[currentnode]]$index)<=neps){
      tree[[currentnode]]$left=NA
      tree[[currentnode]]$right=NA
      
      #其他情况需要计算lossfunc(y|x)
    }else{
      
      #查找使得系数最小的最优特征和最优切分点
      loss=apply(x[tree[[currentnode]]$index,],2,lossfunc,y[tree[[currentnode]]$index])
      variable=which.min(loss[1,])
      
      #当没有更多可供选择的特征时，loss(y|x)=Inf
      #所以若min(loss(y|x))=Inf则代表没有更多可供选择的特征，此时设置当前节点为叶节点
      if (loss[1,variable]==Inf){
        tree[[currentnode]]$left=NA
        tree[[currentnode]]$right=NA
        
        #其他情况当前节点均非叶节点
      }else{
        
        #更新最优特征和最优切分点
        tree[[currentnode]]$variable=variable
        tree[[currentnode]]$divchoice=loss[2,variable]
        
        #更新当前节点的左右子节点，新建子节点
        #左子节点为满足分支条件的数据,右子节点为不满足分支条件的数据
        tree[[currentnode]]$left=length(tree)+1
        tree[[length(tree)+1]]=treenode(tree[[currentnode]]$index[x[tree[[currentnode]]$index,tree[[currentnode]]$variable]<=tree[[currentnode]]$divchoice],
                                        No=length(tree)+1,father=currentnode)
        tree[[currentnode]]$right=length(tree)+1
        tree[[length(tree)+1]]=treenode(tree[[currentnode]]$index[x[tree[[currentnode]]$index,tree[[currentnode]]$variable]>tree[[currentnode]]$divchoice],
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
  if(!cut) return(tree) else return(cart.r.cut(x,y,tree,alphaeps))
}

#CART回归树剪枝算法
cart.r.cut=function(x,y,tree,alphaeps){
  alpha=rep(NA,length(tree))
  #对每一个子树计算满足L(father)=L(tree)的alpha
  for (i in 1:length(tree)){
    if (!is.na(tree[[i]]$left)){
      alpha[i]=sum((y[tree[[i]]$index]-tree[[i]]$value)^2)-sum((y[tree[[tree[[i]]$left]]$index]-tree[[tree[[i]]$left]]$value)^2)-sum((y[tree[[tree[[i]]$right]]$index]-tree[[tree[[i]]$right]]$value)^2)
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
cart.r.predict=function(xpre,tree){
  
  #对单条记录进行预测
  cart.r.predictONE=function(xpre,tree){
    #从根节点开始查找CART分类树
    currentnode=1
    #当其存在子节点时，若满足分支条件，则将当前节点更新左子节点，否则更新为右子节点
    while(!is.na(tree[[currentnode]]$left)){
      if (xpre[tree[[currentnode]]$variable]<=tree[[currentnode]]$divchoice){
        currentnode=tree[[currentnode]]$left
      }else{currentnode=tree[[currentnode]]$right}
    }
    return(tree[[currentnode]]$value)
  }
  
  #对xpre用apply
  if (nrow(xpre)==1){return(cart.r.predictONE(xpre,tree))}else{
    return(apply(xpre,1,cart.r.predictONE,tree))
  }
}

#test part
x=matrix(c(1:10,rep(1:5,2)),nrow = 10)
x=as.data.frame(x)
y=c(4.5,4.75,4.91,5.34,5.8,7.05,7.9,8.23,8.7,9)

tree=cart.r(x,y,neps=2)
cart.r.cut(x,y,tree,alphaeps=0.2)
treecut=cart.r(x,y,neps=2,cut=T,alphaeps=0.2)

xpre=x[2,]
cart.r.predict(xpre,tree)
cart.r.predict(x,tree)
cart.r.predict(xpre,treecut)
cart.r.predict(x,treecut)
