#0.加载包

library(class)
library(stringr)
library(dplyr)
#library(progress)
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
  Train_test_frame$Label<-factor(Train_test_frame$Label,
                                 levels = c(0,1))
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
  }}

#开始测试合适的模型
#建立空白的结果表格
OB_Model_result_frame<-data.frame(matrix(nrow = 0,
                                         ncol = 10))
colnames(OB_Model_result_frame)<-c('Train_rows','Probs','Test_labels','K',
                                   'AUC','number_1','number_0',
                                   'Cor','MAE','ENPEP_Result')
#提取ENPEP数据
ENPEP_data<-exp_change[which(exp_change$Symbol=='ENPEP'),]
rownames(ENPEP_data)<-'ENPEP'
ENPEP_data<-ENPEP_data[,-1]
ENPEP_data<-ENPEP_data[-which(colnames(ENPEP_data)=='Label')]

#循环测试的次数
All_cycle_times<-10000
jishu<-0
#pb <- progress_bar$new(total =
#                         All_cycle_times*(OB_train_num+
#                                         Fat_train_num))
#建立大循环：随机选取固定数量训练集和测试机
for (cycle_times in 1:All_cycle_times) {
  temp_cycle_times<-cycle_times
  #随机选取训练基因
  OB_rows<-sample(which(Train_test_frame$Label==1),
                  OB_train_num,
                  replace = FALSE)
  N_rows<-sample(which(Train_test_frame$Label==0),
                 Fat_train_num,
                 replace = FALSE)
  train_frame<-Train_test_frame[c(OB_rows,N_rows),]
  ob_train_label<-train_frame$Label
  train_frame<-train_frame[,-which(colnames(train_frame)=='Label')]
  #随机选取测试基因
  test_frame<-Train_test_frame[-c(OB_rows,N_rows),]
  ob_test_label<-test_frame$Label
  test_frame<-test_frame[,-which(colnames(test_frame)=='Label')]
  #建立模型
  pred <- knn(train = train_frame,
              test = test_frame,
              cl = ob_train_label,
              k = 116,prob = T)
  ROC_frame<-data.frame(pre=as.numeric(attr(pred,"prob")),
                        answer=as.numeric(as.vector(ob_test_label)))
  p1<-roc(answer ~ pre, ROC_frame,
          level =c(0,1),direction='<')
  #检测ENPEP结果
  ENPEP_Result<-knn(train = train_frame,
                    test = ENPEP_data,
                    cl = ob_train_label,
                    k = 116,prob = T)
  #结果写入表格
  jishu<-jishu+1
  OB_Model_result_frame[jishu,]<-c(numbers_to_char(c(OB_rows,N_rows)),
                                   numbers_to_char(as.numeric(attr(pred,"prob"))),
                                   numbers_to_char(as.vector(ob_test_label)),
                                   116,
                                   as.numeric(p1$auc),
                                   length(which(ob_train_label==1)),
                                   length(which(ob_train_label==0)),
                                   cor(as.numeric(attr(pred,"prob")),
                                       as.numeric(as.vector(ob_test_label))),
                                   MAE(as.numeric(as.vector(ob_test_label)),
                                       as.numeric(attr(pred,"prob"))),
                                   as.vector(ENPEP_Result))
  #pb$tick()  # 更新进度条
}

#D.输出训练结果
write.csv(OB_Model_result_frame,'OB_Model_result_frame_KNN_K116.csv')


