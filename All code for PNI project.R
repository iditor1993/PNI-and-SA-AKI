# 调用当前工作目录下的脚本(阈值效应函数)
source("23final.R")
library(tidyverse)
library(readxl)
raw_data <- read_xlsx("mimiciv_export.xlsx")
complete_data <- raw_data %>% select(-c(subject_id,hadm_id,stay_id,marital_status,insurance,
                                        language,dischtime,admittime,icu_intime,icu_outtime,
                                        careunit,dod,hospital_expire_flag,icu_expire_flag,survival_days,
                                        smoker,alcohol_abuse,mild_liver_disease,peripheral_vascular_disease,dementia,
                                        peptic_ulcer_disease,paraplegia,metastatic_solid_tumor,aids,diagnosis_long_title1,
                                        meld,sepsis,sepsis_time,aki,aki_time,lvef_min,microbiology,
                                        sbp_ni_first,dbp_ni_first,mbp_ni_first,ph_first,po2_first,pco2_first,
                                        so2_first,aado2_first,aado2_calc_first,baseexcess_first,totalco2_first,
                                        methemoglobin_first,carboxyhemoglobin_first,calcium_free_first,
                                        nrbc_first,eosinophils_abs_first,basophils_abs_first,eosinophils_first,
                                        monocytes_first,basophils_first,monocytes_abs_first,lymphocytes_first,
                                        metamyelocytes_first,immature_granulocytes_first,atypical_lymphocytes_first,
                                        hematocrit_first,mch_first,mcv_first,mchc_first,rdw_first,rdwsd_first,rbc_first,
                                        crp_first,total_protein_first,globulin_first,fibrinogen_first,thrombin_first,
                                        d_dimer_first,ck_mb_first,troponin_t_first,ntprobnp_first,alp_first,ggt_first,
                                        ck_cpk_first,amylase_first,ld_ldh_first,glucose_first,hemoglobin_a1c_first,
                                        HDL_first,LDL_first,antibiotic_first,cam_icu_first,heart_rhythm_first,
                                        invasiveVent_hours_first,non_invasiveVent_hours_first,HFNC_hours_first,
                                        input_rbc_sum,input_plt_sum,input_plasma_sum))
complete_data <- as.data.frame(complete_data)
#将性别变量进行转换0=Female，1=Male
complete_data <- complete_data %>% 
  mutate(
    gender = factor(ifelse(gender == "F", 0, 1),
                        levels = 0:1,
                        labels = c("female","male"))
  )
#将年龄、住院总时间、ICU住院时间、身高转换为数值型并取整
complete_data$admission_age <- round(as.numeric(complete_data$admission_age))
complete_data$los_hospital <- round(as.numeric(complete_data$los_hospital))
complete_data$los_icu <- round(as.numeric(complete_data$los_icu))
complete_data$height <- round(as.numeric(complete_data$height))

#将体温、淋巴细胞绝对值、中性粒细胞绝对值、CVP值转换为数值型
complete_data$temperature_first <- as.numeric(complete_data$temperature_first)
complete_data$lymphocytes_abs_first <- as.numeric(complete_data$lymphocytes_abs_first)
complete_data$neutrophils_abs_first <- as.numeric(complete_data$neutrophils_abs_first)
complete_data$CVP_first <- as.numeric(complete_data$CVP_first)

#对人种数据进行合并
library(tidyverse)
library(forcats)
complete_data <- complete_data %>% 
  mutate(
    race = str_squish(str_to_upper(race)),     # 去空格、统一大写
    race = fct_recode(race,
                      "White"    = "WHITE",
                      "White"    = "WHITE - RUSSIAN",
                      "White"    = "WHITE - EASTERN EUROPEAN",
                      "White"    = "WHITE - OTHER EUROPEAN",
                      "White"    = "WHITE - BRAZILIAN",
                      "Black"    = "BLACK/AFRICAN AMERICAN",
                      "Black"    = "BLACK/CARIBBEAN ISLAND",
                      "Black"    = "BLACK/AFRICAN",
                      "Black"    = "BLACK/CAPE VERDEAN",
                      "Asian"    = "ASIAN",
                      "Asian"    = "ASIAN - CHINESE",
                      "Asian"    = "ASIAN - KOREAN",
                      "Asian"    = "ASIAN - ASIAN INDIAN",
                      "Asian"    = "ASIAN - SOUTH EAST ASIAN",
                      "Hispanic" = "HISPANIC/LATINO - CENTRAL AMERICAN",
                      "Hispanic" = "HISPANIC/LATINO - COLUMBIAN",
                      "Hispanic" = "HISPANIC/LATINO - CUBAN",
                      "Hispanic" = "HISPANIC/LATINO - DOMINICAN",
                      "Hispanic" = "HISPANIC/LATINO - GUATEMALAN",
                      "Hispanic" = "HISPANIC/LATINO - HONDURAN",
                      "Hispanic" = "HISPANIC/LATINO - MEXICAN",
                      "Hispanic" = "HISPANIC/LATINO - PUERTO RICAN",
                      "Hispanic" = "HISPANIC/LATINO - SALVADORAN",
                      "Hispanic" = "HISPANIC OR LATINO",
                      "Other"    = "OTHER"
    ),
    race = fct_other(race, keep = c("White","Black","Asian","Hispanic")),  # 剩余→Other
    race = factor(race, levels = c("White","Black","Asian","Hispanic","Other"))
  )

# 查看结果
table(complete_data$race)


#剔除入院时缺失白蛋白与淋巴细胞计数的患者
library(dplyr)
complete_data <- complete_data %>%
  filter(!is.na(albumin_first),     # 剔除 albumin_first 缺失的样本
         !is.na(lymphocytes_abs_first))  # 剔除 lymphocytes_abs_first 缺失的样本

#计算每个样本入院时的预后营养指数（PNI）
library(dplyr)
complete_data <- complete_data %>%
  mutate(PNI = 10 * albumin_first + 5 * lymphocytes_abs_first)

#根据传统的 3 × IQR（四分位距）规则（位于 [Q1–3 × IQR，Q3 + 3 × IQR] 之外的值）来识别和消除异常值
## 1. 计算 PNI 的四分位距与传统 3×IQR 边界
pni_quartiles <- quantile(complete_data$PNI, probs = c(0.25, 0.75), na.rm = TRUE)
Q1 <- pni_quartiles[1]
Q3 <- pni_quartiles[2]
IQR_val <- Q3 - Q1

lower_bound <- Q1 - 3 * IQR_val
upper_bound <- Q3 + 3 * IQR_val

cat(sprintf("PNI 3×IQR 边界: [%.2f, %.2f]\n", lower_bound, upper_bound))

## 2. 生成逻辑向量标记异常值（TRUE=异常）
outlier_flag <- complete_data$PNI < lower_bound | complete_data$PNI > upper_bound

## 3. 剔除异常值，得到新数据框
complete_data_clean <- complete_data[!outlier_flag, ]

## 4. 查看清理结果
cat(sprintf("原始样本量: %d\n", nrow(complete_data)))
cat(sprintf("剔除异常值数量: %d\n", sum(outlier_flag, na.rm = TRUE)))
cat(sprintf("剩余样本量: %d\n", nrow(complete_data_clean)))


#缺失值观察与可视化
library(naniar)
missing_summary_mimic <- miss_var_summary(complete_data_clean)#显示各列缺失值占比
print(missing_summary_mimic)
gg_miss_var(complete_data_clean) + theme_bw()#缺失值可视化

#保留缺失值≤20%的变量
library(dplyr)
# 1. 计算每个变量的缺失比例
miss_rate <- sapply(complete_data_clean, function(x) mean(is.na(x)))
# 2. 找出缺失率 <= 20% 的变量名
keep_var <- names(miss_rate)[miss_rate <= 0.2]
# 3. 只保留这些变量
complete_data_clean <- complete_data_clean %>% select(all_of(keep_var))
# 可选：查看被剔除的变量
cat("剔除的变量（缺失率>20%）：",
    paste(setdiff(names(miss_rate), keep_var), collapse = ", "), "\n")

#对表格数据缺失值进行多重插补
library(mice)#使用mice包进行数据多重插补
set.seed(123)
mimic_com <- complete_data_clean
mimic_com <- mice(mimic_com,m=5,maxit = 10,method = "pmm",seed = 123,print = TRUE)#设置种子，输出过程
mimic_com <- mice::complete(mimic_com)

#对数据集的变量名进行重命名
mimic_com<- mimic_com %>%
  rename_with(~sub("_first$", "", .x))   # 去掉所有变量名末尾的 _first
mimic_com <- mimic_com %>% 
  rename(
    age = admission_age,
    immunity_inhibitor = custom_IMMUINHIB,
    glucocorticoid = custom_GC,
    albumin_infusion = custom_Albtaken,
    AKI_stage = aki_score,
    MV = InvasiveVent,
    GCS = gcs_dynamic,
    SOFA = sofa_dynamic,
    APS_III = apsiii,
    SAPS_II = sapsii,
    SIRS = sirs
  )

#对使用血管活性药物的变量进行转换
library(dplyr)
mimic_com <- mimic_com %>%
  mutate(
    vasoactive_agent = ifelse((custom_VP + custom_EPIs + custom_DOPA + custom_DOBU) == 0, 0, 1)
  )

#对PNI的四分位数进行分组
Q <- quantile(mimic_com$PNI, probs = c(.25, .5, .75), na.rm = TRUE)
mimic_com$group <- with(mimic_com,
                        1 + (PNI > Q[1]) +
                          (PNI > Q[2]) +
                          (PNI > Q[3]))

# 查看分组结果
table(mimic_com$group)

#使用tableone包进行基线值比较并生成表格
# 0. 安装并加载所需包（首次运行需安装）
# install.packages(c("tableone", "dplyr", "writexl"))
library(tableone)
library(dplyr)

# 1. 设定分组变量与需要汇总的变量
#    这里把 group 设为分组，其余变量全部放进汇总
#    如果你的数据里还有非分析变量（如 ID），请从 vars 里剔除
drop_vars <- c("group","custom_VP","custom_EPIs","custom_DOPA","custom_DOBU")
vars  <- setdiff(names(mimic_com), drop_vars)   # 自动去掉分组变量及其他不需要的变量
group <- "group"                              # 分组变量名

# 2. 告诉 tableone 哪些变量是连续型、哪些是分类型
#    连续型默认用正态分布描述（mean±sd），如想改用 median[IQR]：
#    把对应变量名放进 nonnormal 向量即可
nonnormal  <- c("age","los_hospital","los_icu","weight","APS_III",
                "SAPS_II","SIRS","lods","charlson","oasis","GCS",
                "SOFA","temperature","heart_rate","resp_rate",
                "sbp","dbp","mbp","spo2","wbc","lymphocytes_abs",
                "neutrophils_abs","neutrophils","platelet","albumin",
                "creatinine","bun","aniongap","calcium_total","pt",
                "ptt","inr","alt","ast","bilirubin_total","chloride",
                "potassium","sodium","bicarbonate","hemoglobin","magnesium",
                "phosphate","input_all_sum","output_all_sum","output_urine_sum",
                "PNI")   # 举例：这两个用 median[IQR]

cat_vars   <- c("gender","race","mortality_28","mortality_30","mortality_90",
                "mortality_365","hypertension" ,"diabetes","myocardial_infarct",
                "congestive_heart_failure","severe_liver_disease","renal_disease",
                "chronic_pulmonary_disease","cerebrovascular_disease","malignant_cancer",
                "rheumatic_disease","AKI_stage","MV","cpr","crrt","immunity_inhibitor",
                "glucocorticoid","albumin_infusion","vasoactive_agent")  # 举例：这两个是分类型（0/1 或因子）

## 1. 先单独算一份 Overall（不分层）
tab_overall <- CreateTableOne(
  vars      = vars,
  data      = mimic_com,
  factorVars= cat_vars,
  test      = FALSE   # Overall 不需要 P 值
)

## 2. 再算分层 + P 值（4 组）
tab_strata <- CreateTableOne(
  vars      = vars,
  strata    = group,
  data      = mimic_com,
  factorVars= cat_vars,
  test      = TRUE
)

## 3. 分别导出为矩阵
mat_overall <- print(tab_overall,
                     nonnormal   = nonnormal, # 指定用 median[IQR] 的变量
                     showAllLevels = TRUE,    # 分类变量显示全水平
                     printToggle = FALSE,     # 不打印，只返回矩阵
                     noSpaces    = TRUE,      # 控制台对齐
                     varLabels = TRUE         # 如果有变量标签会一起显示
                     )      
mat_strata  <- print(tab_strata,
                     nonnormal   = nonnormal,
                     showAllLevels = TRUE,
                     printToggle = FALSE,
                     noSpaces    = TRUE)



## 5. 分别写 CSV
write.csv(mat_overall, file = "Table1_with_overall.csv", na = "")
write.csv(mat_strata, file = "Table1_with_strata.csv", na = "")


#对PNI分组患者生成30天死亡KM图
library(survival)
library(survminer)
#生成生存曲线用数据与表格
mimic_com_km <- mimic_com %>% select(los_hospital,mortality_30,mortality_90,mortality_365,group)
mimic_com_km <- mimic_com_km %>% mutate(hospstay_30d = ifelse(los_hospital>= 30,30,los_hospital))
diff <- survdiff(Surv(hospstay_30d,mortality_30)~group,data = mimic_com_km)
pvalue <- 1-pchisq(diff$chisq,df=(4-1))#df自由度
if(pvalue<0.001){
  pvalue="Log-rank test \n p<0.001"
}else{
  pvalue=paste0("Log-rank test","\n","p=",sprintf("%.03f",pvalue))
}
fit <- survfit(Surv(hospstay_30d,mortality_30)~group,data = mimic_com_km)
#绘制
mimic_com_30d <- ggsurvplot(fit,
                             data = mimic_com_km,
                             conf.int = F,
                             #pval = T,
                             pval = pvalue,
                             pval.size=5,
                             pval.coord = c(0, 0.5),  # 调整P值位置，x=0, y=0.5
                             legend.labs=c("Quartile 1","Quartile 2","Quartile 3","Quartile 4"),
                             #legend.labs=levels(factor(mimic_km[,"group"])),
                             legend.title="PNI Group",
                             xlab="Time(day)",
                             break.time.by=5,
                             risk.table.title="Number at risk",
                             risk.table = T,
                             risk.table.height=.25,
                             censor = T,
                             censor.shape=124,
                             censor.size=2,
                             ggtheme = theme_bw(),
                             ylim = c(0.4, 1),
                             type = "spline",  # 使生存曲线更圆滑
                             palette = c("red", "blue","darkgreen","orange")  # 设置曲线颜色
)  
print(mimic_com_30d)
jpeg("Figure survival_plot_30d.jpeg", width = 8, height = 6, units = "in", res = 600)
print(mimic_com_30d)
dev.off()

#绘制四个组别住院时间箱线图
library(tidyverse)
library(ggpubr)
mimic_com_boxplot <- mimic_com %>% select(group,mortality_30,los_hospital,los_icu)
mimic_com_boxplot <- mimic_com_boxplot %>% mutate(los_hospital_90d = ifelse(los_hospital>=95,95,los_hospital))
mimic_com_boxplot$group <- factor(mimic_com_boxplot$group,
                                  levels = c("1","2","3","4"),
                                  labels = c("Quartile 1","Quartile 2","Quartile 3","Quartile 4"))#将分组设置标签

mimic_hosp_d <- ggboxplot(mimic_com_boxplot,x="group",y="los_hospital_90d",color = "group",
                          palette = c("red","blue","darkgreen","orange"),
                          ggtheme = theme_bw(),
                          legend.title="Group",
                          legend.labs=c("Quartile 1","Quartile 2","Quartile 3","Quartile 4"),
                          xlab = "",
                          ylab = "Time(day)",
                          ylim = c(0,100)
                          #main= "MIMIC-IV 数据库中谵妄组和非谵妄组的 ICU 住院时间箱线图"
)+
  #stat_compare_means()#添加Wilcoxon秩和检验计算的P值+
  #stat_compare_means(label = "p.signif",
  #symnum.args = list(cutpoints=c(0,0.0001,0.001,0.01,0.05,Inf),
  #symbols=c("<0.0001","<0.001","<0.01","<0.05","ns")))
  stat_compare_means(
    method = "kruskal.test",   # 4 组整体比较
    label = "p.signif",
    symnum.args = list(cutpoints=c(0,0.01,0.05,Inf),
                      symbols=c("P<0.01","P<0.05","ns")),#获取区间并设置区间的P值怎么显示
    label.x = "Quartile 1",#确定P值的x轴位置
    label.y = 96,         #适当抬高，避免压住箱线
    vjust = 0             # 0 = 底对齐，1 = 顶对齐
    )
print(mimic_hosp_d)
jpeg("Figure boxplot.jpeg", width = 8, height = 6, units = "in", res = 600)
print(mimic_hosp_d)
dev.off()

#进行COX分析
# 1. 进行单因素cox分析
library(survival)
library(dplyr)
library(tidyverse)
library(broom)
library(stringr)

#对数据集进行整理，去除不需要的变量
mimic_com_cox <- mimic_com %>% select(-c(los_icu,mortality_28,mortality_90,mortality_365,
                                         custom_VP,custom_EPIs,custom_DOPA,custom_DOBU))
mimic_com_cox$group <- as.factor(mimic_com_cox$group) 
mimic_com_cox$AKI_stage <- as.factor(mimic_com_cox$AKI_stage)

# 加载所需的包
library(survival)
library(tidyverse)

# 执行单因素Cox回归分析
univ_cox_models <- function(data, time_var, status_var) {
  # 获取所有候选变量（排除生存时间和状态变量）
  candidate_vars <- names(data)[!names(data) %in% c(time_var, status_var)]
  
  results <- list()
  
  for (var in candidate_vars) {
    # 构建公式
    formula <- as.formula(paste("Surv(", time_var, ", ", status_var, ") ~ ", var))
    
    # 拟合Cox模型
    cox_model <- coxph(formula, data = data)
    
    # 提取模型摘要
    model_summary <- summary(cox_model)
    
    # 获取结果
    if (nrow(model_summary$coefficients) == 1) {
      # 连续变量或二分类变量
      coef <- model_summary$coefficients[1, "coef"]
      hr <- model_summary$coefficients[1, "exp(coef)"]
      hr_ci_lower <- model_summary$conf.int[1, "lower .95"]
      hr_ci_upper <- model_summary$conf.int[1, "upper .95"]
      p_value <- model_summary$coefficients[1, "Pr(>|z|)"]
      
      results[[var]] <- data.frame(
        Variable = var,
        Levels = NA,
        Coef = coef,
        HR = hr,
        HR_CI_Lower = hr_ci_lower,
        HR_CI_Upper = hr_ci_upper,
        P_Value = p_value,
        stringsAsFactors = FALSE
      )
    } else {
      # 多分类变量（多个水平）
      for (i in 1:nrow(model_summary$coefficients)) {
        level_name <- rownames(model_summary$coefficients)[i]
        coef <- model_summary$coefficients[i, "coef"]
        hr <- model_summary$coefficients[i, "exp(coef)"]
        hr_ci_lower <- model_summary$conf.int[i, "lower .95"]
        hr_ci_upper <- model_summary$conf.int[i, "upper .95"]
        p_value <- model_summary$coefficients[i, "Pr(>|z|)"]
        
        results[[paste(var, level_name, sep = "_")]] <- data.frame(
          Variable = var,
          Levels = level_name,
          Coef = coef,
          HR = hr,
          HR_CI_Lower = hr_ci_lower,
          HR_CI_Upper = hr_ci_upper,
          P_Value = p_value,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # 合并所有结果
  do.call(rbind, results)
}

# 执行分析
cox_results <- univ_cox_models(mimic_com_cox, "los_hospital", "mortality_30")

# 查看结果
print(cox_results)

# 保存为CSV文件
write.csv(cox_results, "univariable_cox_analysis_results_01.csv", row.names = FALSE)

#查看PNI分组的整体趋势，进行单因素COX回归
mimic_com_cox$group <- as.numeric(mimic_com_cox$group)#将因子型变量转换为数值型变量，统计中用来观察整体趋势
#进行单因素cox回归
cox_model_group <- coxph(Surv(los_hospital, mortality_30) ~ group, data = mimic_com_cox)
summary(cox_model_group)#得出的P值就是p for trend
mimic_com_cox$group <- as.factor(mimic_com_cox$group)#将变量转换回因子型

# 可选：美化P值显示（科学计数法转小数）
cox_results_pretty <- cox_results %>%
  mutate(P_Value = format.pval(P_Value, digits = 3, eps = 0.001),
         HR_95CI = paste0(round(HR_CI_Lower, 3), "-", round(HR_CI_Upper, 3)),
         HR = round(HR, 3)) %>%
  select(Variable, Levels, Coef, HR, HR_95CI, P_Value)

# 保存美化后的结果
write.csv(cox_results_pretty, "univariable_cox_analysis_results_pretty.csv", row.names = FALSE)

#保留单因素COX回归P＜0.05的变量
sig_df <- cox_results %>% filter(P_Value < 0.05)      # 保留整行
sig_var <- cox_results %>% filter(P_Value < 0.05) %>% pull(Variable)  # 只要 Variable 向量

#对单因素COX回归有意义的结果进行多因素COX回归
library(dplyr)
mimic_com_multi_cox <- mimic_com_cox %>% 
  select(los_hospital,mortality_30,all_of(sig_var))   # 确保 sig_var 是字符向量，用 all_of 防止拼写警告

# 精简版本 - 包含系数
library(survival)

# 多因素COX回归
covariates <- setdiff(names(mimic_com_multi_cox), c("los_hospital", "mortality_30"))
cox_model_multi <- coxph(Surv(los_hospital, mortality_30) ~ ., data = mimic_com_multi_cox)

# 提取结果 - 包含系数
results <- data.frame(
  variable = names(cox_model_multi$coefficients),
  coefficient = cox_model_multi$coefficients,  # 原始系数
  HR = exp(cox_model_multi$coefficients),
  CI_lower = exp(confint(cox_model_multi)[, 1]),
  CI_upper = exp(confint(cox_model_multi)[, 2]),
  p_value = summary(cox_model_multi)$coefficients[, 5]
)

# 保存结果
write.csv(results, "cox_results_with_coefficients.csv", row.names = FALSE)





#纳入相关协变量进行调整观察PNI与死亡的关联
#model 2 纳入基线合并症与干预措施观察PNI与患者30天死亡的关联
names(mimic_com_cox)#获取变量名
mimic_cox_model2 <- mimic_com_cox %>% select(los_hospital,mortality_30,age,weight,myocardial_infarct,congestive_heart_failure,
                                             severe_liver_disease,malignant_cancer,cpr,crrt,PNI)
#PNI连续值在model2中的HR值
cox_multi_model2 <- coxph(Surv(los_hospital, mortality_30) ~ ., data = mimic_cox_model2)

# 提取结果 - 包含系数
results_cox_model2 <- data.frame(
  variable = names(cox_multi_model2$coefficients),
  coefficient = cox_multi_model2$coefficients,  # 原始系数
  HR = exp(cox_multi_model2$coefficients),
  CI_lower = exp(confint(cox_multi_model2)[, 1]),
  CI_upper = exp(confint(cox_multi_model2)[, 2]),
  p_value = summary(cox_multi_model2)$coefficients[, 5]
)

# 保存结果
write.csv(results_cox_model2, "cox_multi_model2_with_coefficients_PNIconti.csv", row.names = FALSE)


#新数据包含group，无连续PNI
mimic_cox_model2 <- mimic_com_cox %>% select(los_hospital,mortality_30,age,weight,myocardial_infarct,congestive_heart_failure,
                                             severe_liver_disease,malignant_cancer,cpr,crrt,group)

#对group变量转换为数值型，观察PNI分组整体趋势
mimic_cox_model2$group <- as.numeric(mimic_cox_model2$group)


#对新数据集进行多因素COX回归分析(p for trend)
cox_multi_model2 <- coxph(Surv(los_hospital, mortality_30) ~ ., data = mimic_cox_model2)

# 提取结果 - 包含系数
results_cox_model2 <- data.frame(
  variable = names(cox_multi_model2$coefficients),
  coefficient = cox_multi_model2$coefficients,  # 原始系数
  HR = exp(cox_multi_model2$coefficients),
  CI_lower = exp(confint(cox_multi_model2)[, 1]),
  CI_upper = exp(confint(cox_multi_model2)[, 2]),
  p_value = summary(cox_multi_model2)$coefficients[, 5]
)

# 保存结果
write.csv(results_cox_model2, "cox_multi_model2_with_coefficients_ptrend.csv", row.names = FALSE)

#对group变量转换为因子型，观察PNI分组的具体变化
mimic_cox_model2$group <- as.factor(mimic_cox_model2$group)
#对新数据集进行多因素COX回归分析(PNI分组)
cox_multi_model2 <- coxph(Surv(los_hospital, mortality_30) ~ ., data = mimic_cox_model2)

# 提取结果 - 包含系数
results_cox_model2 <- data.frame(
  variable = names(cox_multi_model2$coefficients),
  coefficient = cox_multi_model2$coefficients,  # 原始系数
  HR = exp(cox_multi_model2$coefficients),
  CI_lower = exp(confint(cox_multi_model2)[, 1]),
  CI_upper = exp(confint(cox_multi_model2)[, 2]),
  p_value = summary(cox_multi_model2)$coefficients[, 5]
)

# 保存结果
write.csv(results_cox_model2, "cox_multi_model2_with_coefficients.csv", row.names = FALSE)


#纳入相关协变量进行调整观察PNI与死亡的关联
#model 3 纳入基线合并症、干预措施、实验室参数与疾病严重评分观察PNI与患者30天死亡的关联
names(mimic_com_cox)#获取变量名
mimic_cox_model3 <- mimic_com_cox %>% select(los_hospital,mortality_30,age,weight,myocardial_infarct,congestive_heart_failure,
                                             severe_liver_disease,malignant_cancer,cpr,crrt,platelet,hemoglobin,
                                             creatinine,bun,aniongap,calcium_total,pt,ptt,inr,ast,bilirubin_total,
                                             potassium,bicarbonate,magnesium,phosphate,SOFA,GCS,charlson,AKI_stage,PNI)

#PNI连续值在model3中的HR值
cox_multi_model3 <- coxph(Surv(los_hospital, mortality_30) ~ ., data = mimic_cox_model3)

# 提取结果 - 包含系数
results_cox_model3 <- data.frame(
  variable = names(cox_multi_model3$coefficients),
  coefficient = cox_multi_model3$coefficients,  # 原始系数
  HR = exp(cox_multi_model3$coefficients),
  CI_lower = exp(confint(cox_multi_model3)[, 1]),
  CI_upper = exp(confint(cox_multi_model3)[, 2]),
  p_value = summary(cox_multi_model3)$coefficients[, 5]
)

# 保存结果
write.csv(results_cox_model3, "cox_multi_model3_with_coefficients_PNIconti.csv", row.names = FALSE)

#新数据包含group，无连续PNI
mimic_cox_model3 <- mimic_com_cox %>% select(los_hospital,mortality_30,age,weight,myocardial_infarct,congestive_heart_failure,
                                             severe_liver_disease,malignant_cancer,cpr,crrt,platelet,hemoglobin,
                                             creatinine,bun,aniongap,calcium_total,pt,ptt,inr,ast,bilirubin_total,
                                             potassium,bicarbonate,magnesium,phosphate,SOFA,GCS,charlson,AKI_stage,group)


#对group变量转换为数值型，观察PNI分组整体趋势
mimic_cox_model3$group <- as.numeric(mimic_cox_model3$group)

#对新数据集进行多因素COX回归分析(p for trend)
cox_multi_model3 <- coxph(Surv(los_hospital, mortality_30) ~ ., data = mimic_cox_model3)

# 提取结果 - 包含系数
results_cox_model3 <- data.frame(
  variable = names(cox_multi_model3$coefficients),
  coefficient = cox_multi_model3$coefficients,  # 原始系数
  HR = exp(cox_multi_model3$coefficients),
  CI_lower = exp(confint(cox_multi_model3)[, 1]),
  CI_upper = exp(confint(cox_multi_model3)[, 2]),
  p_value = summary(cox_multi_model3)$coefficients[, 5]
)

# 保存结果
write.csv(results_cox_model3, "cox_multi_model3_with_coefficients_ptrend.csv", row.names = FALSE)


#对group变量转换为因子型，观察PNI分组的具体变化
mimic_cox_model3$group <- as.factor(mimic_cox_model3$group)
#对新数据集进行多因素COX回归分析(PNI分组)
cox_multi_model3 <- coxph(Surv(los_hospital, mortality_30) ~ ., data = mimic_cox_model3)

# 提取结果 - 包含系数
results_cox_model3 <- data.frame(
  variable = names(cox_multi_model3$coefficients),
  coefficient = cox_multi_model3$coefficients,  # 原始系数
  HR = exp(cox_multi_model3$coefficients),
  CI_lower = exp(confint(cox_multi_model3)[, 1]),
  CI_upper = exp(confint(cox_multi_model3)[, 2]),
  p_value = summary(cox_multi_model3)$coefficients[, 5]
)

# 保存结果
write.csv(results_cox_model3, "cox_multi_model3_with_coefficients.csv", row.names = FALSE)

#RCS回归分析及图
library(tidyverse)
library(ggrcs)
library(rms)
library(ggplot2)
library(scales)
library(cowplot)

# 设置rms包所需的数据分布
mimic_rcs <- mimic_com_cox %>% select(los_hospital,mortality_30,age,weight,myocardial_infarct,congestive_heart_failure,
                                      severe_liver_disease,malignant_cancer,cpr,crrt,platelet,hemoglobin,
                                      creatinine,bun,aniongap,calcium_total,pt,ptt,inr,ast,bilirubin_total,
                                      potassium,bicarbonate,magnesium,phosphate,SOFA,charlson,AKI_stage,
                                      PNI,group)
PNI_dist <- datadist(mimic_rcs)
options(datadist ="PNI_dist")

#拟合Cox比例风险回归模型
fit <- cph(Surv(los_hospital, mortality_30 == 1) ~ rcs(PNI, 3) + age + weight + myocardial_infarct + 
             congestive_heart_failure + severe_liver_disease + malignant_cancer + cpr + crrt + 
             platelet + hemoglobin + creatinine + bun + aniongap + calcium_total + pt + ptt + 
             inr + ast + bilirubin_total + potassium + bicarbonate + magnesium + phosphate + 
             SOFA + charlson + AKI_stage, x = TRUE, y = TRUE, data = mimic_rcs)

# 提取P值
anova_result <- anova(fit)
p_overall <- format.pval(anova_result["PNI","P"], digits = 3, eps = 0.001)
p_nonlinear <- format.pval(anova_result[" Nonlinear","P"], digits = 3, eps = 0.001)
p_text <- paste("P for overall =", p_overall,"\nP for nonlinear =", p_nonlinear)

# 绘制单一RCS曲线
rcs_plot1 <- ggrcs(
  data = mimic_rcs,
  fit = fit,
  x ="PNI",
  histcol ="#1F77B4", # JAMA风格深蓝色
  histbinwidth = 3,
  ribcol ="#BC3C29", # JAMA风格浅灰色
  ribalpha = 0.3,
  xlab ="PNI",
  ylab ="Hazard Ratio of 30-day mortality",
  title ="Association Between PNI and 30-day Mortality Hazard Ratio") +
  #annotate("text", x = min(mimic_rcs$PNI), y = Inf, label = p_text, # 手动添加P值，左上角
  #         hjust = 0, vjust = 1.2, size = 4, color ="black") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face ="bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color ="black")
  )


# 显示绘图
print(rcs_plot1)

jpeg("Figure RCS.jpeg", width = 10, height = 8, units = "in", res = 600)
print(rcs_plot1)
dev.off()


#寻找PNI的阈值
library(survival)
# 注意：这里拟合的模型不含RCS项
fit_coxph <- cph(Surv(los_hospital, mortality_30 == 1) ~ PNI + age + weight + myocardial_infarct + 
                     congestive_heart_failure + severe_liver_disease + malignant_cancer + cpr + crrt + 
                     platelet + hemoglobin + creatinine + bun + aniongap + calcium_total + pt + ptt + 
                     inr + ast + bilirubin_total + potassium + bicarbonate + magnesium + phosphate + 
                     SOFA + charlson + AKI_stage, data = mimic_rcs) # 补齐其他协变量
# 假设 cut.tab 函数已定义或通过源文件加载
threshold_result <- cuttab(fit_coxph, "PNI", mimic_rcs)
print(threshold_result)
write.csv(threshold_result,"threshold result table.csv",row.names = FALSE)

#在图中添加阈值参考线（可选）
rcs_plot1 + geom_vline(aes(xintercept=34.85),colour="#BB0000",linetype="dashed")

jpeg("Figure RCS_final.jpeg", width = 10, height = 8, units = "in", res = 600)
rcs_plot1 + geom_vline(aes(xintercept=34.85),colour="#BB0000",linetype="dashed")
dev.off()

#制作交互作用森林图
#加载R包
#install.packages("devtools")
library(devtools)
#devtools::install_github("adayim/forestploter")
#install.packages("jstable")
library(jstable)
library(survival)
library(grid)
library(forestploter)
library(tidyverse)
library(labelled)
#针对PNI的阈值进行分组
mimic_cox_p_intera <- mimic_com_cox %>% 
  select(los_hospital,mortality_30,age,hypertension,diabetes,myocardial_infarct,
        congestive_heart_failure,renal_disease,chronic_pulmonary_disease,
        cerebrovascular_disease,rheumatic_disease,severe_liver_disease,
        malignant_cancer,cpr,crrt,AKI_stage,PNI)

mimic_cox_p_intera <- mimic_cox_p_intera %>% mutate(PNI_group=ifelse(PNI < 34.85,0,1))
mimic_cox_p_intera <- mimic_cox_p_intera %>% mutate(age_group=ifelse(age < 65,0,1))
#mimic_cox_p_intera <- mimic_cox_p_intera %>% filter(PNI_group == 0)

mimic_cox_p_intera <- mimic_cox_p_intera %>% 
  mutate(age_group=factor(age_group,levels = c(0,1),labels = c("<65 years","≥65 years")),
         hypertension=factor(hypertension,levels = c(0,1),labels = c("No","Yes")),
         diabetes=factor(diabetes,levels = c(0,1),labels = c("No","Yes")),
         myocardial_infarct=factor(myocardial_infarct,levels = c(0,1),labels = c("No","Yes")),
         congestive_heart_failure=factor(congestive_heart_failure,levels = c(0,1),labels = c("No","Yes")),
         renal_disease=factor(renal_disease,levels = c(0,1),labels = c("No","Yes")),
         chronic_pulmonary_disease=factor(chronic_pulmonary_disease,levels = c(0,1),labels = c("No","Yes")),
         cerebrovascular_disease=factor(cerebrovascular_disease,levels = c(0,1),labels = c("No","Yes")),
         rheumatic_disease=factor(rheumatic_disease,levels = c(0,1),labels = c("No","Yes")),
         severe_liver_disease=factor(severe_liver_disease,levels = c(0,1),labels = c("No","Yes")),
         malignant_cancer=factor(malignant_cancer,levels = c(0,1),labels = c("No","Yes")),
         cpr=factor(cpr,levels = c(0,1),labels = c("No","Yes")),
         crrt=factor(crrt,levels = c(0,1),labels = c("No","Yes")),
         AKI_stage=factor(AKI_stage,levels = c(1,2,3),labels = c("Stage 1","Stage 2","Stage 3")),
         PNI_group=factor(PNI_group,levels = c(0,1),labels = c("Low PNI","High PNI"))
  )

mimic_cox_p_intera <- mimic_cox_p_intera %>% 
  set_variable_labels(age_group ="Age",
                      hypertension = "Hypertension",
                      diabetes = "Diabetes",
                      renal_disease = "Renal disease",
                      myocardial_infarct = "Myocardial infarct",
                      congestive_heart_failure = "Congestive heart failure",
                      chronic_pulmonary_disease = "Chronic pulmonary disease",
                      cerebrovascular_disease = "Cerebrovascular disease",
                      rheumatic_disease = "Rheumatic disease",
                      severe_liver_disease = "Severe liver disease",
                      malignant_cancer = "Malignant cancer",
                      cpr = "CPR",
                      crrt = "CRRT",
                      AKI_stage = "AKI stage")


res_PNI <- TableSubgroupMultiCox(formula = Surv(los_hospital, mortality_30) ~ PNI_group, #提供Cox模型的公式（本示例中我们对PNI值对于生存的影响进行Cox回归分析）
                                 data = mimic_cox_p_intera,
                                 var_subgroups = c("age_group","hypertension","diabetes","renal_disease",
                                                   "myocardial_infarct","congestive_heart_failure",
                                                   "chronic_pulmonary_disease","cerebrovascular_disease",
                                                   "rheumatic_disease","severe_liver_disease",
                                                   "malignant_cancer","cpr","crrt","AKI_stage"), #指定需要进行亚组分析的变量名称
                                 count_by = "PNI_group",  #按照不同的PNI值对人群计数
                                 line = F,    #是否在不同亚组变量间添加一条新的线段
                                 event = T    #按照不同的治疗方式对人群计数
)

res_PNI <- res_PNI %>%  
  mutate(across(everything(), ~ifelse(is.na(.), " ", .)),
         ` ` = paste(rep(" ", 20), collapse = " "),
         `Hazard Ratio (95% CI)` = ifelse(`Point Estimate` == " ", " ", paste0(`Point Estimate`, " (", Lower, "-", Upper, ")", sep = "")),
         across(c(`Point Estimate`,Lower,Upper), ~as.numeric(.))) |> 
  rename("Subgroup" = "Variable",
         "Low PNI" = "Count(PNI_group=Low PNI)",
         "High PNI" = "Count(PNI_group=High PNI)")

res_PNI$Count <- gsub("\\(.*?\\)", "", res_PNI$Count)
res_PNI$'Low PNI' <- gsub("\\(.*?\\)", "", res_PNI$'Low PNI')
res_PNI$'High PNI' <- gsub("\\(.*?\\)", "", res_PNI$'High PNI')

view(res_PNI)

# Define theme
tm <- forest_theme(
  base_size = 10,        # 设置文本的基础大小
  
  colhead=list(fg_params=list(hjust=0)), # 设置列名对齐方式 
  core = list(padding = unit(c(4, 4), "mm")), #增加森林图每行的宽度
  
  # 设置置信区间的外观
  ci_pch = 16,           # 置信区间点的形状
  ci_col = "black",    # 置信区间的边框颜色
  ci_fill = "black",      # 置信区间的填充颜色
  ci_alpha = 0.8,        # 置信区间的透明度
  ci_lty = 1,            # 置信区间的线型
  ci_lwd = 1.0,          # 置信区间的线宽
  ci_Theight = 0.2,      # 设置T字在置信区间末端的高度，默认是NULL
  
  # 设置参考线的外观
  refline_gp = gpar(lwd = 1,            # 参考线的线宽
                    lty = "dashed",     # 参考线的线型
                    col = "grey20"      # 参考线的颜色
  ),
  
  # 设置垂直线的外观
  vertline_lwd = 1,         # 垂直线的线宽，可以添加一条额外的垂直线，如果没有就不显示
  vertline_lty = "dashed",  # 垂直线的线型
  vertline_col = "grey20",  # 垂直线的颜色
  
  # 设置脚注的字体大小、字体样式和颜色
  footnote_gp = gpar(
    # fontfamily = "italic",   # 脚注文本的字体
    cex = 0.6,               # 脚注字体大小
    fontface = "plain",      # 脚注字体样式
    col = "black"            # 脚注文本的颜色
  ),
)

#res_PNI_no_overall <- res_PNI[-1,]

p <- forest(res_PNI[,c(1, 3, 4, 13, 14, 12)],
            est = res_PNI$`Point Estimate`,
            lower = res_PNI$Lower, 
            upper = res_PNI$Upper,
            sizes = res_PNI$`Point Estimate`,
            ci_column = 4,
            ref_line = 1,
            arrow_lab = c("High PNI Better", "Low PNI Better"),
            xlim = c(0, 2),
            ticks_at = c(0,0.5, 1, 1.5, 2),
            #footnote = "This is the demo data. Please feel free to change\nanything you want.",
            theme = tm)
p <- insert_text(p,
                 text = "no. of patients with events/total no.",
                 row = 1, col = 2:3,
                 just = "left",
                 gp = gpar(cex = 0.7, col = "black", fontface = "italic"))
p

jpeg("Figure forestplot.jpeg", width = 8, height = 13, units = "in", res = 600)
print(p)
dev.off()

#注意，主效应为连续变量的最好还是做COX回归交互作用的交互效应图，暂时没有找到做的很好的办法，留着以后做

#计算C指数
# 加载必要的包
library(survival)
library(rms)
library(pROC)
library(ggplot2)
library(tidyverse)
#整理统计用数据集
mimic_c_index <- mimic_com_cox %>% select(los_hospital,mortality_30,APS_III,SAPS_II,
                                          SIRS,lods,oasis,SOFA,
                                          charlson,PNI,group)
mimic_c_index <- mimic_c_index %>% rename(PNI_group = group)

#计算C指数
# 建立计算C-index函数
calculate_c_index <- function(predictor, outcome) {
  roc_obj <- roc(outcome, predictor, quiet = TRUE)
  return(auc(roc_obj))
}

# 1. 30天死亡率预测 - 连续PNI
cat("\n=== 30天死亡率预测 (连续PNI) ===\n")

c_index_30day <- data.frame(
  Model = character(),
  C_index = numeric(),
  stringsAsFactors = FALSE
)

# 单独模型
models_30day <- list(
  "PNI" = mimic_c_index$PNI,
  "OASIS" = mimic_c_index$oasis,
  "SOFA" = mimic_c_index$SOFA,
  "APSIII" = mimic_c_index$APS_III,
  "SAPSII" = mimic_c_index$SAPS_II,
  "SIRS"= mimic_c_index$SIRS,
  "LODS"= mimic_c_index$lods,
  "Charlson"=mimic_c_index$charlson
)

for (model_name in names(models_30day)) {
  c_index <- calculate_c_index(models_30day[[model_name]], mimic_c_index$mortality_30)
  c_index_30day <- rbind(c_index_30day, 
                         data.frame(Model = model_name, C_index = c_index))
}

# 组合模型
combined_models <- list(
  "PNI + OASIS" = glm(mortality_30 ~ PNI + oasis, data = mimic_c_index, family = binomial)$fitted,
  "PNI + SOFA" = glm(mortality_30 ~ PNI + SOFA, data = mimic_c_index, family = binomial)$fitted,
  "PNI + APSIII" = glm(mortality_30 ~ PNI + APS_III, data = mimic_c_index, family = binomial)$fitted,
  "PNI + SAPSII" = glm(mortality_30 ~ PNI + SAPS_II, data = mimic_c_index, family = binomial)$fitted,
  "PNI + SIRS" = glm(mortality_30 ~ PNI + SIRS, data = mimic_c_index, family = binomial)$fitted,
  "PNI + LODS" = glm(mortality_30 ~ PNI + lods, data = mimic_c_index, family = binomial)$fitted,
  "PNI + Charlson" = glm(mortality_30 ~ PNI + charlson, data = mimic_c_index, family = binomial)$fitted
)

for (model_name in names(combined_models)) {
  c_index <- calculate_c_index(combined_models[[model_name]], mimic_c_index$mortality_30)
  c_index_30day <- rbind(c_index_30day, 
                         data.frame(Model = model_name, C_index = c_index))
}

print(c_index_30day)
write.csv(c_index_30day,"c_index_30d_PNI_continue.csv",row.names = FALSE)


# 2. 30天死亡率预测 - PNI四分位数
cat("\n=== 30天死亡率预测 (PNI四分位数) ===\n")

# 将分类变量转换为数值评分
mimic_c_index$PNI_group <- as.numeric(mimic_c_index$PNI_group)


c_index_30day_Q <- data.frame(
  Model = character(),
  C_index = numeric(),
  stringsAsFactors = FALSE
)

# 单独模型
models_30day_Q <- list(
  "PNI_group" = mimic_c_index$PNI_group,
  "OASIS" = mimic_c_index$oasis,
  "SOFA" = mimic_c_index$SOFA,
  "APSIII" = mimic_c_index$APS_III,
  "SAPSII" = mimic_c_index$SAPS_II,
  "SIRS"= mimic_c_index$SIRS,
  "LODS"= mimic_c_index$lods,
  "Charlson"=mimic_c_index$charlson
)

for (model_name in names(models_30day_Q)) {
  c_index_Q <- calculate_c_index(models_30day_Q[[model_name]], mimic_c_index$mortality_30)
  c_index_30day_Q <- rbind(c_index_30day_Q, 
                           data.frame(Model = model_name, C_index = c_index_Q))
}

# 组合模型
combined_models_Q <- list(
  "PNI group + OASIS" = glm(mortality_30 ~ PNI_group + oasis, data = mimic_c_index, family = binomial)$fitted,
  "PNI group + SOFA" = glm(mortality_30 ~ PNI_group + SOFA, data = mimic_c_index, family = binomial)$fitted,
  "PNI group + APSIII" = glm(mortality_30 ~ PNI_group + APS_III, data = mimic_c_index, family = binomial)$fitted,
  "PNI group + SAPSII" = glm(mortality_30 ~ PNI_group + SAPS_II, data = mimic_c_index, family = binomial)$fitted,
  "PNI group + SIRS" = glm(mortality_30 ~ PNI_group + SIRS, data = mimic_c_index, family = binomial)$fitted,
  "PNI group + LODS" = glm(mortality_30 ~ PNI_group + lods, data = mimic_c_index, family = binomial)$fitted,
  "PNI group + Charlson" = glm(mortality_30 ~ PNI_group + charlson, data = mimic_c_index, family = binomial)$fitted
)

for (model_name in names(combined_models_Q)) {
  c_index_Q <- calculate_c_index(combined_models_Q[[model_name]], mimic_c_index$mortality_30)
  c_index_30day_Q <- rbind(c_index_30day_Q, 
                           data.frame(Model = model_name, C_index = c_index_Q))
}

print(c_index_30day_Q)
write.csv(c_index_30day_Q,"c_index_30d_PNI_quantile.csv",row.names = FALSE)


# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)

# 读取数据
c_index_data <- read.csv("c_index_30d_PNI_continue.csv")

# 创建分组数据框
plot_data <- data.frame(
  Group = rep(c("OASIS", "SOFA", "APSIII", "SAPSII", "SIRS", "LODS", "Charlson"), each = 3),
  Subgroup = rep(c("PNI Alone", "Score Alone", "Score+PNI"), 7),
  C_index = c(
    # OASIS组
    c_index_data$C_index[c_index_data$Model == "PNI"],
    c_index_data$C_index[c_index_data$Model == "OASIS"],
    c_index_data$C_index[c_index_data$Model == "PNI + OASIS"],
    
    # SOFA组
    c_index_data$C_index[c_index_data$Model == "PNI"],
    c_index_data$C_index[c_index_data$Model == "SOFA"],
    c_index_data$C_index[c_index_data$Model == "PNI + SOFA"],
    
    # APSIII组
    c_index_data$C_index[c_index_data$Model == "PNI"],
    c_index_data$C_index[c_index_data$Model == "APSIII"],
    c_index_data$C_index[c_index_data$Model == "PNI + APSIII"],
    
    # SAPSII组
    c_index_data$C_index[c_index_data$Model == "PNI"],
    c_index_data$C_index[c_index_data$Model == "SAPSII"],
    c_index_data$C_index[c_index_data$Model == "PNI + SAPSII"],
    
    # SIRS组
    c_index_data$C_index[c_index_data$Model == "PNI"],
    c_index_data$C_index[c_index_data$Model == "SIRS"],
    c_index_data$C_index[c_index_data$Model == "PNI + SIRS"],
    
    # LODS组
    c_index_data$C_index[c_index_data$Model == "PNI"],
    c_index_data$C_index[c_index_data$Model == "LODS"],
    c_index_data$C_index[c_index_data$Model == "PNI + LODS"],
    
    # Charlson组
    c_index_data$C_index[c_index_data$Model == "PNI"],
    c_index_data$C_index[c_index_data$Model == "Charlson"],
    c_index_data$C_index[c_index_data$Model == "PNI + Charlson"]
  )
)

# 设置因子水平以确保正确的顺序
plot_data$Group <- factor(plot_data$Group, 
                          levels = c("OASIS", "SOFA", "APSIII", "SAPSII", "SIRS", "LODS", "Charlson"))
plot_data$Subgroup <- factor(plot_data$Subgroup, 
                             levels = c("PNI Alone", "Score Alone", "Score+PNI"))

# JAMA风格配色
jama_colors <- c("PNI Alone" = "#B24745",    # 红色
                 "Score Alone" = "#00A1D5",   # 蓝色
                 "Score+PNI" = "#80796B")     # 灰色

# 创建柱状图
C_index_plot <- ggplot(plot_data, aes(x = Group, y = C_index, fill = Subgroup)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.9) +
  scale_fill_manual(values = jama_colors, 
                    name = "Model",
                    labels = c("PNI Alone", "Score Alone", "Score+PNI")) +
  labs(#title = "C-index comparison of PNI, clinical scores, and combined models for 30-day mortality prediction",
       x = " ",
       y = "C-index") +
  theme_classic() +  # 使用classic主题，它默认只有左边框和下边框
  theme(
    # 图例位置调整到右上角空白区域
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.9),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    # 图例背景和边框透明
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.key.size = unit(0.5, "cm"),
    
    # 坐标轴和文本样式
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(face = "bold", size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    
    # 标题样式
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, 
                              margin = margin(b = 15)),
    
    # 确保上边框和右边框被移除
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    
    # 边距调整
    plot.margin = margin(20, 20, 20, 20)
  ) +
  # 添加数值标签 - 缩小字号到3.0
  geom_text(aes(label = sprintf("%.3f", C_index)), 
            position = position_dodge(0.8), 
            vjust = -0.5, size = 3.0, fontface = "bold") +
  # 设置y轴范围，确保标签可见
  coord_cartesian(ylim = c(0.55, 0.72)) +
  scale_y_continuous(breaks = seq(0.55, 0.72, by = 0.05))

# 显示图形
print(C_index_plot)

# 保存为高分辨率图片
jpeg("C_index_comparison_JAMA_style.jpeg", width = 12, height = 8, units = "in", res = 600)
print(C_index_plot)
dev.off()
#ggsave("C_index_comparison_JAMA_style.jpeg", plot = C_index_plot, width = 12, height = 8, dpi = 600, bg = "white")



# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)

# 读取数据
c_index_data <- read.csv("c_index_30d_PNI_quantile.csv")

# 创建分组数据框
plot_data_Q <- data.frame(
  Group = rep(c("OASIS", "SOFA", "APSIII", "SAPSII", "SIRS", "LODS", "Charlson"), each = 3),
  Subgroup = rep(c("PNI Quartile Alone", "Score Alone", "Score+PNI Quartile"), 7),
  C_index = c(
    # OASIS组
    c_index_data$C_index[c_index_data$Model == "PNI_group"],
    c_index_data$C_index[c_index_data$Model == "OASIS"],
    c_index_data$C_index[c_index_data$Model == "PNI group + OASIS"],
    
    # SOFA组
    c_index_data$C_index[c_index_data$Model == "PNI_group"],
    c_index_data$C_index[c_index_data$Model == "SOFA"],
    c_index_data$C_index[c_index_data$Model == "PNI group + SOFA"],
    
    # APSIII组
    c_index_data$C_index[c_index_data$Model == "PNI_group"],
    c_index_data$C_index[c_index_data$Model == "APSIII"],
    c_index_data$C_index[c_index_data$Model == "PNI group + APSIII"],
    
    # SAPSII组
    c_index_data$C_index[c_index_data$Model == "PNI_group"],
    c_index_data$C_index[c_index_data$Model == "SAPSII"],
    c_index_data$C_index[c_index_data$Model == "PNI group + SAPSII"],
    
    # SIRS组
    c_index_data$C_index[c_index_data$Model == "PNI_group"],
    c_index_data$C_index[c_index_data$Model == "SIRS"],
    c_index_data$C_index[c_index_data$Model == "PNI group + SIRS"],
    
    # LODS组
    c_index_data$C_index[c_index_data$Model == "PNI_group"],
    c_index_data$C_index[c_index_data$Model == "LODS"],
    c_index_data$C_index[c_index_data$Model == "PNI group + LODS"],
    
    # Charlson组
    c_index_data$C_index[c_index_data$Model == "PNI_group"],
    c_index_data$C_index[c_index_data$Model == "Charlson"],
    c_index_data$C_index[c_index_data$Model == "PNI group + Charlson"]
  )
)

# 设置因子水平以确保正确的顺序
plot_data_Q$Group <- factor(plot_data_Q$Group, 
                            levels = c("OASIS", "SOFA", "APSIII", "SAPSII", "SIRS", "LODS", "Charlson"))
plot_data_Q$Subgroup <- factor(plot_data_Q$Subgroup, 
                               levels = c("PNI Quartile Alone", "Score Alone", "Score+PNI Quartile"))

# JAMA风格配色
jama_colors <- c("PNI Quartile Alone" = "#B24745",    # 红色
                 "Score Alone" = "#00A1D5",   # 蓝色
                 "Score+PNI Quartile" = "#80796B")     # 灰色

# 创建柱状图
C_index_quartile_plot <- ggplot(plot_data_Q, aes(x = Group, y = C_index, fill = Subgroup)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.9) +
  scale_fill_manual(values = jama_colors, 
                    name = "Model",
                    labels = c("PNI Quartile Alone", "Score Alone", "Score+PNI Quartile")) +
  labs(#title = "C-index comparison of PNI quartile, clinical scores, and combined models for 30-day mortality prediction",
       x = " ",
       y = "C-index") +
  theme_classic() +  # 使用classic主题，它默认只有左边框和下边框
  theme(
    # 图例位置调整到右上角空白区域
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.90),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    # 图例背景和边框透明
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.key.size = unit(0.5, "cm"),
    
    # 坐标轴和文本样式
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(face = "bold", size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    
    # 标题样式
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, 
                              margin = margin(b = 15)),
    
    # 确保上边框和右边框被移除
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    
    # 边距调整
    plot.margin = margin(20, 20, 20, 20)
  ) +
  # 添加数值标签 - 缩小字号到3.0
  geom_text(aes(label = sprintf("%.3f", C_index)), 
            position = position_dodge(0.8), 
            vjust = -0.5, size = 3.0, fontface = "bold") +
  # 设置y轴范围，确保标签可见
  coord_cartesian(ylim = c(0.55, 0.72)) +
  scale_y_continuous(breaks = seq(0.55, 0.72, by = 0.05))

# 显示图形
print(C_index_quartile_plot)

# 保存为高分辨率图片
jpeg("C_index_quartile_comparison_JAMA_style.jpeg", width = 12, height = 8, units = "in", res = 600)
print(C_index_quartile_plot)
dev.off()
#ggsave("C_index_quartile_comparison_JAMA_style.jpeg", plot = C_index_quantile_plot, width = 12, height = 8, dpi = 600, bg = "white")
