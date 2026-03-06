
library(shiny)
library(survival)

# =========================
# 读取模型与数据
# =========================
final_model <- readRDS("final_model.rds")
train_data  <- readRDS("train_data.rds")

final_vars <- c("hscrp", "zhongxing", "putaotang", "qiudanbai",
                "ALT", "CEA", "rs2279199", "BMI",
                "general", "differentiated", "combination", "surgery", "Rs11615")

# =========================
# 中文标签
# =========================
var_labels <- c(
  hscrp         = "超敏C反应蛋白（hs-CRP）",
  zhongxing     = "中性粒细胞计数",
  putaotang     = "葡萄糖",
  qiudanbai     = "球蛋白",
  ALT           = "丙氨酸氨基转移酶（ALT）",
  CEA           = "癌胚抗原（CEA）",
  rs2279199     = "UMPS（rs2279199）",
  BMI           = "BMI分组",
  general       = "肿瘤大体分型",
  differentiated= "肿瘤分化程度",
  combination   = "是否联用靶向药物",
  surgery       = "肿瘤手术史",
  Rs11615       = "ERCC1（rs11615）"
)

# 连续变量与分类变量分组
continuous_vars <- c("hscrp", "zhongxing", "putaotang", "qiudanbai", "ALT", "CEA")
categorical_vars <- c("rs2279199", "BMI", "general", "differentiated", "combination", "surgery", "Rs11615")

# =========================
# 样式
# =========================
app_css <- "
body {
  background-color: #f5f7fb;
  font-family: 'Microsoft YaHei', 'PingFang SC', sans-serif;
}
.main-header {
  background: linear-gradient(90deg, #1f4e79, #2f75b5);
  color: white;
  padding: 18px 24px;
  border-radius: 12px;
  margin-bottom: 20px;
}
.main-header h2 {
  margin: 0;
  font-weight: 700;
}
.main-header p {
  margin: 6px 0 0 0;
  opacity: 0.95;
}
.panel-box {
  background: white;
  border-radius: 14px;
  padding: 18px;
  box-shadow: 0 4px 12px rgba(0,0,0,0.08);
  margin-bottom: 18px;
}
.section-title {
  font-size: 18px;
  font-weight: 700;
  color: #1f3c5b;
  margin-bottom: 12px;
  border-left: 5px solid #2f75b5;
  padding-left: 10px;
}
.predict-btn {
  background-color: #2f75b5 !important;
  border-color: #2f75b5 !important;
  color: white !important;
  font-weight: bold;
  width: 100%;
  height: 48px;
  border-radius: 10px;
}
.card-row {
  display: flex;
  gap: 16px;
  flex-wrap: wrap;
}
.value-card {
  flex: 1;
  min-width: 180px;
  background: linear-gradient(180deg, #ffffff, #f7fbff);
  border-radius: 14px;
  padding: 18px;
  box-shadow: 0 3px 10px rgba(0,0,0,0.07);
  border-top: 5px solid #2f75b5;
}
.value-card h4 {
  margin: 0;
  font-size: 16px;
  color: #355c7d;
}
.value-card .big-number {
  font-size: 28px;
  font-weight: 700;
  color: #1f3c5b;
  margin-top: 10px;
}
.risk-box-low {
  background: #edf7ed;
  color: #256029;
  border-left: 6px solid #4caf50;
  padding: 14px;
  border-radius: 10px;
  font-size: 16px;
  font-weight: 600;
}
.risk-box-high {
  background: #fff3f3;
  color: #a12622;
  border-left: 6px solid #e53935;
  padding: 14px;
  border-radius: 10px;
  font-size: 16px;
  font-weight: 600;
}
.small-note {
  color: #6b7280;
  font-size: 13px;
}
"

# =========================
# 工具函数
# =========================

# 自动生成输入控件
make_input_ui <- function(var, data) {
  x <- data[[var]]
  label_cn <- ifelse(!is.na(var_labels[var]), var_labels[var], var)
  
  if (is.factor(x)) {
    selectInput(
      inputId = var,
      label = label_cn,
      choices = levels(x),
      selected = levels(x)[1]
    )
  } else if (is.character(x)) {
    choices <- unique(x)
    choices <- choices[!is.na(choices)]
    selectInput(
      inputId = var,
      label = label_cn,
      choices = choices,
      selected = choices[1]
    )
  } else if (is.logical(x)) {
    selectInput(
      inputId = var,
      label = label_cn,
      choices = c("FALSE", "TRUE"),
      selected = "FALSE"
    )
  } else {
    med <- stats::median(x, na.rm = TRUE)
    step_value <- ifelse(all(x %% 1 == 0, na.rm = TRUE), 1, 0.1)
    
    numericInput(
      inputId = var,
      label = label_cn,
      value = round(med, 2),
      step = step_value
    )
  }
}

# 构建 newdata
build_newdata <- function(input, vars, train_data) {
  res <- vector("list", length(vars))
  names(res) <- vars
  
  for (v in vars) {
    x <- train_data[[v]]
    
    if (is.factor(x)) {
      res[[v]] <- factor(input[[v]], levels = levels(x))
    } else if (is.character(x)) {
      res[[v]] <- as.character(input[[v]])
    } else if (is.logical(x)) {
      res[[v]] <- as.logical(input[[v]])
    } else if (is.integer(x)) {
      res[[v]] <- as.integer(input[[v]])
    } else {
      res[[v]] <- as.numeric(input[[v]])
    }
  }
  
  as.data.frame(res, stringsAsFactors = FALSE)
}

# 计算训练集风险分层阈值：以线性预测值中位数为cut-off
train_lp <- predict(final_model, newdata = train_data, type = "lp")
risk_cutoffs <- quantile(train_lp, probs = c(0.33, 0.66), na.rm = TRUE)

cutoff_low  <- risk_cutoffs[1]
cutoff_high <- risk_cutoffs[2]

# 提取指定时间点生存率
extract_surv_prob <- function(fit, time_point) {
  s <- summary(fit, times = time_point)
  if (length(s$surv) == 0) return(NA_real_)
  as.numeric(s$surv[1])
}

# 卡片HTML
make_card <- function(title, value_text) {
  div(
    class = "value-card",
    h4(title),
    div(class = "big-number", value_text)
  )
}

# =========================
# UI
# =========================
ui <- fluidPage(
  tags$head(
    tags$style(HTML(app_css))
  ),
  
  div(
    class = "main-header",
    h2("结直肠癌患者无进展生存个体化临床预测系统"),
    p("基于Cox回归模型的无进展生存期（PFS）风险评估与个体化预测")
  ),
  
  fluidRow(
    column(
      width = 4,
      
      div(
        class = "panel-box",
        div(class = "section-title", "患者信息输入"),
        p(class = "small-note", "请根据患者实际情况填写连续变量和分类变量后点击“开始预测”。"),
        
        div(class = "section-title", style = "font-size:16px;", "连续变量"),
        uiOutput("continuous_inputs"),
        
        tags$hr(),
        
        div(class = "section-title", style = "font-size:16px;", "分类变量"),
        uiOutput("categorical_inputs"),
        
        br(),
        actionButton("predict_btn", "开始预测", class = "predict-btn")
      )
    ),
    
    column(
      width = 8,
      
      div(
        class = "panel-box",
        div(class = "section-title", "预测结果概览"),
        htmlOutput("risk_group_ui"),
        br(),
        uiOutput("summary_cards"),
        br(),
        tableOutput("prob_table")
      ),
      
      div(
        class = "panel-box",
        div(class = "section-title", "个体化无进展生存生存曲线"),
        plotOutput("surv_plot", height = "420px")
      ),
      
      div(
        class = "panel-box",
        div(class = "section-title", "模型输出指标"),
        verbatimTextOutput("risk_text")
      )
    )
  )
)

# =========================
# Server
# =========================
server <- function(input, output, session) {
  
  # 连续变量输入
  output$continuous_inputs <- renderUI({
    ui_list <- lapply(continuous_vars, make_input_ui, data = train_data)
    do.call(tagList, ui_list)
  })
  
  # 分类变量输入
  output$categorical_inputs <- renderUI({
    ui_list <- lapply(categorical_vars, make_input_ui, data = train_data)
    do.call(tagList, ui_list)
  })
  
  # 构造患者数据
  patient_data <- eventReactive(input$predict_btn, {
    build_newdata(input, final_vars, train_data)
  })
  
  # 模型结果
  pred_result <- eventReactive(input$predict_btn, {
    newdata <- patient_data()
    
    lp <- predict(final_model, newdata = newdata, type = "lp")
    risk <- predict(final_model, newdata = newdata, type = "risk")
    fit <- survfit(final_model, newdata = newdata)
    
    pfs_9  <- extract_surv_prob(fit, 9)
    pfs_12 <- extract_surv_prob(fit, 12)
    pfs_15 <- extract_surv_prob(fit, 15)
    
    risk_group <- ifelse(
      lp <= cutoff_low,
      "低风险",
      ifelse(lp <= cutoff_high, "中风险", "高风险")
    )
    
    list(
      lp = lp,
      risk = risk,
      fit = fit,
      pfs_9 = pfs_9,
      pfs_12 = pfs_12,
      pfs_15 = pfs_15,
      risk_group = risk_group
    )
  })
  
  # 风险分层提示
  output$risk_group_ui <- renderUI({
    res <- pred_result()
    
    if (res$risk_group == "高风险") {
      div(class = "risk-box-high",
          paste0("风险分层结果：", res$risk_group,
                 "。提示疾病进展风险较高，建议加强随访监测。"))
      
    } else if (res$risk_group == "中风险") {
      div(style="background:#fff7e6;border-left:6px solid #ffa940;padding:14px;border-radius:10px;font-weight:600;",
          paste0("风险分层结果：", res$risk_group,
                 "。提示中等进展风险，建议常规随访评估。"))
      
    } else {
      div(class = "risk-box-low",
          paste0("风险分层结果：", res$risk_group,
                 "。提示进展风险较低，但仍需定期随访。"))
    }
  })
  
  # 9/12/15月预测卡片
  output$summary_cards <- renderUI({
    res <- pred_result()
    
    p9  <- ifelse(is.na(res$pfs_9),  "NA", paste0(round(res$pfs_9 * 100, 1), "%"))
    p12 <- ifelse(is.na(res$pfs_12), "NA", paste0(round(res$pfs_12 * 100, 1), "%"))
    p15 <- ifelse(is.na(res$pfs_15), "NA", paste0(round(res$pfs_15 * 100, 1), "%"))
    
    div(
      class = "card-row",
      make_card("9个月无进展生存概率",  p9),
      make_card("12个月无进展生存概率", p12),
      make_card("15个月无进展生存概率", p15)
    )
  })
  
  # 表格结果
  output$prob_table <- renderTable({
    res <- pred_result()
    
    data.frame(
      时间点 = c("9个月", "12个月", "15个月"),
      `预测PFS概率` = c(
        ifelse(is.na(res$pfs_9),  NA, sprintf("%.4f", res$pfs_9)),
        ifelse(is.na(res$pfs_12), NA, sprintf("%.4f", res$pfs_12)),
        ifelse(is.na(res$pfs_15), NA, sprintf("%.4f", res$pfs_15))
      )
    )
  }, striped = TRUE, bordered = TRUE, spacing = "m")
  
  # 生存曲线
  output$surv_plot <- renderPlot({
    res <- pred_result()
    
    plot(
      res$fit,
      xlab = "随访时间（月）",
      ylab = "无进展生存概率",
      main = "个体化无进展生存预测曲线",
      lwd = 3,
      mark.time = TRUE,
      conf.int = TRUE,
      col = "#2f75b5",
      xaxs = "i",
      yaxs = "i"
    )
    
    abline(v = c(9, 12, 15), lty = 2, col = "gray60")
    grid()
  })
  
  # 模型数值输出
  output$risk_text <- renderPrint({
    res <- pred_result()
    
    cat("线性预测值（LP）:", round(res$lp, 4), "\n")
    cat("相对风险（Risk Score）:", round(res$risk, 4), "\n")
    cat("风险分层阈值：\n")
    cat("低风险阈值 (33%):", round(cutoff_low,4), "\n")
    cat("高风险阈值 (66%):", round(cutoff_high,4), "\n")
    cat("风险分层结果:", res$risk_group, "\n")
  })
}

shinyApp(ui = ui, server = server)
