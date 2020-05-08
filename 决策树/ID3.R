#用于计算信息增益（互信息）的函数
I=function(x,y,base=2){
  n=length(y)
  py=table(y)/n
  Hy=-sum(py*log(py,base))
  
  px=table(x)/n
  pxy=table(x,y)
  pxy=pxy/rowSums(pxy)
  t=log(pxy,base)
  if (any(pxy==0)) t[pxy==0]=0
  Hyx=sum(-rowSums(pxy*t)*px)
  
  return(Hy-Hyx)
}

#定义ID3树节点形式(ID3决策树不一定为二叉树，每一节点的子节点个数等于所选特征xi的取值数)
treenode=function(x,No=NA,chosen=NULL,father=NA,son=NULL,variable=NA,divchoice=NA,value=NA){
  return(list(No=No,chosen=chosen,father=father,son=son,variable=variable,divchoice=divchoice,value=value,index=x))
}

#ID3生成算法
#选择信息增益最大的特征
#叶节点条件：
#1.当前节点的的x所对应的y均属于同一类Ck
#2.特征集为空集
#3.节点的最大信息增益小于阈值(eps)
id3=function(x,y,eps=1e-6,base=2,cut=FALSE,alpha=NULL){
  n=nrow(x)
  p=ncol(x)
  tree=list()
  currentnode=1
  tree[[currentnode]]=treenode(1:n,1)
  
  while (TRUE){
    
    #更新当前节点所对应的类标记
    tree[[currentnode]]$value=names(which.max(table(y[tree[[currentnode]]$index])))
    
    #若当前节点的的x所对应的y均属于同一类Ck或特征集为空集，则将该节点修改为叶节点
    if (all(y[tree[[currentnode]]$index]==y[tree[[currentnode]]$index][1])||
        setequal(tree[[currentnode]]$chosen,1:p)){
      tree[[currentnode]]$son=NA
	  
    #其他情况需要计算信息增益
    }else{
      #查找使得信息增益最大的变量
      if (is.null(tree[[currentnode]]$chosen)){
        IAg=apply(x[tree[[currentnode]]$index,],2,I,y[tree[[currentnode]]$index],base)
      }else{
        IAg=apply(x[tree[[currentnode]]$index,-tree[[currentnode]]$chosen],2,I,y[tree[[currentnode]]$index],base)
      }
      #tree[[tree[[currentnode]]$son]]$variable=which.max(IAg)
      variable=which.max(IAg)
      #将选出的变量添加到chosen中
      tree[[currentnode]]$chosen=c(tree[[currentnode]]$chosen,variable)
      
      #若该节点的最大信息增益小于阈值，则将该节点作为叶节点
      if (IAg[variable]<eps) {
        tree[[currentnode]]$son=0
		
      #其他情况当前节点均非叶节点
      }else{
        #更新当前节点的子节点
        subx=x[tree[[currentnode]]$index,variable]
        k=length(levels(subx))
        tree[[currentnode]]$son=seq(length(tree)+1,by=1,length.out = k)
        #建立当前节点的子节点
        for (i in 1:k){
          tree[[length(tree)+1]]=treenode(tree[[currentnode]]$index[subx==levels(subx)[i]],
                                          No=length(tree)+1,chosen=tree[[currentnode]]$chosen,
                                          father=currentnode,variable=variable,divchoice=levels(subx)[i])
        }
      }
    }
    
    #遍历整个树，查找下一更新节点
    for (j in 1:length(tree)){
      if (is.na(tree[[j]]$value)) {newnode=tree[[j]]$No;break}
    }
    #当所有节点均被更新完，则退出循环
    if (newnode!=currentnode) {currentnode=newnode} else {break}
  }
  if (cut==FALSE) return(tree) else{
    return(id3.cut(x,y,tree,alpha,base))
  }
}

#ID3剪枝算法
id3.cut=function(x,y,tree,alpha,base=2){
  lossfunction=function(y,index,alpha,base=2){ #ID3损失函数
    py=table(y[index])/length(index)
    py=py*log(py,base)
    py[is.na(py)]=0
    return(-sum(py)*length(y)+alpha)
  }
  dellist=NULL
  for (i in length(tree):1){
    if (is.na(tree[[i]]$son) && (!i%in% dellist)){
      #计算父节点的损失函数值
	  father=tree[[i]]$father
      fatherloss=lossfunction(y,index = tree[[father]]$index,alpha,base) 
	  #计算子节点的损失函数值
      son=tree[[father]]$son
      sonloss=0
      for (j in son){
        sonloss=sonloss+lossfunction(y,index = tree[[j]]$index,alpha,base)
      }
	  #当父节点损失函数值<子节点损失函数值，则将父节点变为叶节点
      if (sonloss>fatherloss){
        dellist=c(dellist,tree[[father]]$son)
        tree[[father]]$son=NA
      }
    }
  }
  if (!is.null(dellist)){
    for (i in dellist){
      tree[[i]]=0
    }
  }
  return(tree)
}

#ID3预测
id3.predict=function(xpre,tree){

  #对单条记录进行预测
  id3.predictONE=function(xpre,tree){
    #从根节点开始查找ID3树
    currentnode=1
	#当其存在子节点时，将当前节点更新为对应分类的子节点
    while(!is.na(tree[[currentnode]]$son)[1]){
      for (i in tree[[currentnode]]$son){
        if (xpre[tree[[i]]$variable]==tree[[i]]$divchoice){
          currentnode=i
          break
        }
      }
    }
    return(tree[[currentnode]]$value)
  }
  
  #对xpre用apply
  if (nrow(xpre)==1){
    return(id3.predictONE(xpre,tree))
  }else return(apply(xpre,1,id3.predictONE,tree))
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

tree=id3(x,y)
id3.cut(x,y,tree,alpha = 0.5)
id3(x,y,cut = TRUE,alpha = 0.5)

xpre=x[4,]
id3.predict(xpre,tree)
id3.predict(x,tree)
tree[[2]]$son=NA
id3.predict(xpre,tree)
id3.predict(x,tree)

