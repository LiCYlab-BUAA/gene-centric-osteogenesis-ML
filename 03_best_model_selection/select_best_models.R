#KNN
{
  data<-read.csv('OB_Model_result_frame_KNN_K116.csv')
  data<-data[which(data$ENPEP_Result==1),]
  library(dplyr)
  data<-arrange(data,desc(AUC))
  
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
  roc_knn<-roc(as.numeric(chafen(data$Test_labels[1]))~as.numeric(chafen(data$Probs[1])),
               levels=c(0,1))
  plot(roc_knn,col='red',grid=TRUE,legacy.axes = TRUE,
       thresholds="best",
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
  result<-jieguo(as.numeric(chafen(data$Probs[1])),0.530)
  T_P<-length(which(as.numeric(chafen(data$Test_labels[1]))==1 & result==1))
  F_P<-length(which(as.numeric(chafen(data$Test_labels[1]))==0 & result==1))
  T_N<-length(which(as.numeric(chafen(data$Test_labels[1]))==0 & result==0))
  F_N<-length(which(as.numeric(chafen(data$Test_labels[1]))==1 & result==0))
  c(T_P,T_N,F_P,F_N)
}

#SVM
{
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
}

#LR
{
  #0.加载包
  
  library(class)
  library(stringr)
  library(dplyr)
  library(pROC)
  
  
  #构建必要的其他函数
  {
    numbers_to_char<-function(x){
      temp_char<-''
      for (i1 in 1:length(x)) {
        temp_char<-paste0(temp_char,', ',x[i1])
      }
      final_result<-str_sub(temp_char,3,nchar(temp_char))
      return(final_result)
    }
    char_to_numbers<-function(x){
      temp_numbers<-c()
      temp_sep<-str_locate_all(x,', ')
      temp_start<-c(1,c(as.numeric(temp_sep[[1]][,2])+1))
      temp_end<-c((as.numeric(temp_sep[[1]][,1])-1),nchar(x))
      for (i2 in 1:length(temp_start)) {
        temp_numbers<-c(
          temp_numbers,
          str_sub(x,temp_start[i2],temp_end[i2])
        )
      }
      temp_numbers<-temp_numbers
      return(temp_numbers)
    }
    judge_anwser<-function(answer,to_check){
      assign('result',list())
      FP<-length(which(to_check==1 & answer==0))
      TP<-length(which(to_check==1 & answer==1))
      FN<-length(which(to_check==0 & answer==1))
      TN<-length(which(to_check==0 & answer==0))
      result[[1]]<-FP/(TN+FP)#FPR
      result[[2]]<-TP/(TP+FN)#TPR
      result[[3]]<-FN/(TP+FN)#FNR
      result[[4]]<-TN/(TN+FP)#TNR
      names(result)<-c('FPR','TPR','FNR','TNR')
      return(result)
    }
  }
  
  #Cubist
  {
    OB_Model_result_frame<-read.csv('OB_Model_result_frame_Cubist.csv')
    
    temp_OB_Model_result_frame<-OB_Model_result_frame[which(
      OB_Model_result_frame$MAE==min(OB_Model_result_frame$MAE)),]
    ROC_frame<-data.frame(Probs=as.numeric(char_to_numbers(temp_OB_Model_result_frame$Probs[1])),
                          Labels=char_to_numbers(temp_OB_Model_result_frame$Test_labels[1]))
    p1<-roc(Labels~Probs,ROC_frame,
            level=c(0,1),direction='<')
    plot(p1,col="red",
         print.auc=TRUE)
    
  }
  
  {
    library(ggplot2)
    library(ggrepel)
    MAE_min<-min(OB_Model_result_frame$MAE)
    Selected_Cor<-OB_Model_result_frame$Cor[which(OB_Model_result_frame$MAE==min(OB_Model_result_frame$MAE))]
    
    ggplot(data=OB_Model_result_frame)+
      geom_point(mapping = aes(x=MAE,y=Cor),size=0.1,color='gray')+
      geom_point(mapping = aes(x=MAE_min,y=Selected_Cor),size=0.5,color='red')+
      theme_bw()+
      theme(axis.title = element_text(face = 'bold'))
  }
}

#RF
{
  
  #读取数据
  model_data<-openxlsx::read.xlsx('OB_Model_result_frame_RF_1.xlsx',
                                  sheet = 'OB_Model_result_frame_RF_1')
  ava_model_data<-model_data[which(model_data$ENPEP_Class_Result==1),]
  
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
  roc_rf<-roc(as.numeric(chafen(ava_model_data$Test_labels[1]))~as.numeric(chafen(ava_model_data$Probs[1])),
              levels=c(0,1))
  plot(roc_rf)
  
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
  result<-jieguo(as.numeric(chafen(ava_model_data$Probs[1])),ava_model_data$Best_threshold[1])
  T_P<-length(which(as.numeric(chafen(ava_model_data$Test_labels[1]))==1 & result==1))
  F_P<-length(which(as.numeric(chafen(ava_model_data$Test_labels[1]))==0 & result==1))
  T_N<-length(which(as.numeric(chafen(ava_model_data$Test_labels[1]))==0 & result==0))
  F_N<-length(which(as.numeric(chafen(ava_model_data$Test_labels[1]))==1 & result==0))
  
  c(T_P,T_N,F_P,F_N)
  
  
  plot(roc_rf,
       col='red',add=TRUE,grid=TRUE,legacy.axes = TRUE,
       main="ROC曲线最佳阈值点",
       thresholds="best", # 基于youden指数选择roc曲线最佳阈值点
       print.auc=T,
       print.thres="best")
  
  #查看表示每个预测变量重要性的得分
  #summary(RF_ob_pred_model)
  importance_otu <- RF_ob_pred_model$importance
  head(importance_otu)
  #或者使用函数 importance()
  importance_otu <- data.frame(importance(RF_ob_pred_model), check.names = FALSE)
  head(importance_otu)
  #作图展示 top30 重要的预测变量
  varImpPlot(RF_ob_pred_model, n.var = min(30, nrow(RF_ob_pred_model$importance)),
             main = 'Top 30 - variable importance')
  
}

plot(roc_knn,col='red',legacy.axes = TRUE)
plot(roc_svm,col='blue',add=TRUE)
plot(p1,col="yellow",
     add=TRUE)
plot(roc_rf,
     col='green',add=TRUE)
legend('bottomright',
       legend = c('k-NN AUC = 0.910','SVM AUC = 0.953',
                  'LR AUC = 0.949','RF AUC = 0.910'),
       col=c('red','blue','yellow','green'),
       lty=1)
