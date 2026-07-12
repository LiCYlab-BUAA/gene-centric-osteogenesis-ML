
#读取数据
model_data<-read.csv('OB_Model_result_frame_SVM.csv')
ava_model_data<-model_data[which(model_data$ENPEP_Result==1),]

chafen<-function(char){
  library(stringr)
  char1<-str_replace_all(char,' ','')
  seps<-as.numeric(str_locate_all(char1,',')[[1]][,1])
  start<-c(1,seps+1)
  end<-c(seps-1,nchar(char1))
  result<-c()
  for (i in 1:length(start)) {
    result<-c(result,
              str_sub(char1,(start[i]),(end[i])))
  }
  return(result)
}

library(pROC)
roc_svm<-roc(as.numeric(chafen(ava_model_data$Test_labels[1]))~as.numeric(chafen(ava_model_data$Probs[1])),
            levels=c(0,1))
plot(roc_svm,col='red',add=TRUE,grid=TRUE,legacy.axes = TRUE,
     main="ROC曲线最佳阈值点",
     thresholds="best", # 基于youden指数选择roc曲线最佳阈值点
     print.auc=T,
     print.thres="best")

jieguo<-function(probs,threshold){
  assign('output',c())
  for (i in 1:length(probs)) {
    if(as.numeric(probs[i])>threshold){
      output<-c(output,1)
    }else if(as.numeric(probs[i])<threshold){
      output<-c(output,0)
    }
  }
  return(output)
}
result<-jieguo(as.numeric(chafen(ava_model_data$Probs[1])),0.5240873)
T_P<-length(which(as.numeric(chafen(ava_model_data$Test_labels[1]))==1 & result==1))
F_P<-length(which(as.numeric(chafen(ava_model_data$Test_labels[1]))==0 & result==1))
T_N<-length(which(as.numeric(chafen(ava_model_data$Test_labels[1]))==0 & result==0))
F_N<-length(which(as.numeric(chafen(ava_model_data$Test_labels[1]))==1 & result==0))

c(T_P,T_N,F_P,F_N)


