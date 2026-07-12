#0.加载包

library(class)
library(stringr)
library(dplyr)
#library(progress)
library(Cubist)
library(pROC)

#准备工作
{
  #读取数据
  exp_change<-read.csv('log2_dataset_frame.csv',header = T)
  colnames(exp_change)[1]<-'Symbol'
  all_genes<-unique(exp_change$Symbol)
  
  #建立标签集
  exp_change$Label<-'Unknown'
  Labels<-read.csv('Label_genes.csv',header = T)
  OB_gene_all<-Labels$Symbol[which(Labels$Label==1)]
  raw_negative_genes<-Labels$Symbol[which(Labels$Label==0)]
  exp_change$Label[which(exp_change$Symbol %in% OB_gene_all)]<-1
  ava_OB_genes<-all_genes[which(all_genes %in% OB_gene_all)]
  ava_selected_negative_genes<-intersect(all_genes,
                                         toupper(raw_negative_genes))
  exp_change$Label[which(exp_change$Symbol %in% ava_selected_negative_genes)]<-0
  
  #可以抽取的训练测试矩阵
  Train_test_frame<-exp_change[which(exp_change$Label!='Unknown'),]
  ava_genes<-Train_test_frame$Symbol
  rownames(Train_test_frame)<-Train_test_frame$Symbol
  Train_test_frame<-Train_test_frame[,-1]
  
  #训练集和测试集的基因数量
  OB_train_num<-round(length(ava_OB_genes)*0.8)#从OB基因抽取的训练集基因数
  OB_test_num<-length(ava_OB_genes)-OB_train_num#从OB基因抽取的测试集基因数
  Fat_train_num<-round(length(ava_selected_negative_genes)*0.8)#从Fat基因抽取的训练集基因数
  Fat_test_num<-length(ava_selected_negative_genes)-Fat_train_num#从Fat基因抽取的测试集基因数
  
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
      temp_numbers<-as.numeric(temp_numbers)
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
    MAE<-function(actual,predicted){
      mean(abs(actual-predicted))
    }
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
  }
  
}

#开始测试合适的模型

#建立空白的结果表格
OB_Model_result_frame<-data.frame(matrix(nrow = 0,
                                         ncol = 9))
colnames(OB_Model_result_frame)<-c('Train_rows','Probs','Test_labels',
                                   'AUC','number_1','number_0',
                                   'Cor','MAE','ENPEP_Result')
#提取ENPEP数据
ENPEP_data<-exp_change[which(exp_change$Symbol=='ENPEP'),]
rownames(ENPEP_data)<-'ENPEP'
ENPEP_data<-ENPEP_data[,-1]

model_data<-read.csv('OB_Model_result_frame_Cubist.csv')

train_rows<-as.numeric(chafen(model_data$Train_rows[1]))
train_frame<-Train_test_frame[train_rows,]
ob_train_label<-train_frame$Label
train_frame<-train_frame[,-which(colnames(train_frame)=='Label')]
#随机选取测试基因
test_frame<-Train_test_frame[-train_rows,]
ob_test_label<-test_frame$Label
test_frame<-test_frame[,-which(colnames(test_frame)=='Label')]
#建立回归模型
ob_pred_model<-cubist(x = train_frame,
                      y = as.numeric(ob_train_label))
#进行预测
ob_test_pred<-predict(ob_pred_model,
                      test_frame)
#效能评价
ROC_frame<-data.frame(pre=as.numeric(ob_test_pred),
                      answer=as.numeric(ob_test_label))
p1<-roc(answer ~ pre, ROC_frame,
        level =c(0,1),direction='<')

plot(p1,col='red',add=TRUE,grid=TRUE,legacy.axes = TRUE,
     main="ROC曲线最佳阈值点",
     thresholds="best", # 基于youden指数选择roc曲线最佳阈值点
     print.auc=T,
     print.thres="best")
#检测ENPEP结果
ENPEP_Result<-predict(ob_pred_model,
                      ENPEP_data)

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
result<-jieguo(as.numeric(ob_test_pred),0.991)
T_P<-length(which(as.numeric(ob_test_label)==1 & result==1))
F_P<-length(which(as.numeric(ob_test_label)==0 & result==1))
T_N<-length(which(as.numeric(ob_test_label)==0 & result==0))
F_N<-length(which(as.numeric(ob_test_label)==1 & result==0))
c(T_P,T_N,F_P,F_N)

library(caret)

varImp(ob_pred_model,data= train_frame)
write.csv(data.frame(varImp(ob_pred_model,data= train_frame)),'varImp.csv')
plot(varImp(ob_pred_model))

dotplot(ob_pred_model)
