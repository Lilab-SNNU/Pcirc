library(shiny)
options(shiny.maxRequestSize = 2*1024*1024^2)
options(show.output.on.console = FALSE)
# Define UI ----
ui <- fluidPage(
  title = "Pcirc Running Page",
  h1('Pcirc Pipline',align='center'),
  br(),
  sidebarLayout(
    
    sidebarPanel(width = 3,
                 h3('INTRODUCTION'),
                 p("This page helps to run",span("Pcirc pipline", style = "color:blue"),
                 p('PCirc is a pipeline to predict plant circular RNA (CircRNA) which based on Python3.'),
                 p('It can identify circRNA from a given RNA-seq data by high-throughput.'),
                 p('You can also get more information in github.'),
                 p('(https://github.com/Lilab-SNNU/Pcirc) '),
                   h4('Usage:'),
                 p('Upload all the files and chick the RUN buttom'))),
    
    fluidRow(column(3,
                    helpText('Running thread, value > 0'),
              numericInput(value = 4,
                           inputId = 'thread',
                           label = 'thread',
                           step = 1,
                           min = 1,
                           width = 100),
              br(),
              helpText('NGS file in fastq format'),
              fileInput(inputId = 'fastqfile',label = 'fastq', multiple = TRUE)),
             column(3,
              helpText('GTF file in fastq format'),
              fileInput(inputId = 'gtffile',label = 'gtf'),
              helpText('genome file in fasta format'),
              fileInput(inputId = 'genomefile',label = 'genome'),
              actionButton(inputId = 'runapp',label = 'RUN')),
              column(2,
                     helpText('parameter:'),
                     textOutput('ele')
           )))
   )
  

# Define server logic ----
server <- function(input, output) {

output$ele <- renderText({
  commandstr <- paste('python3 PCirc_shiny.py -g',input$genomefile[[4]],
                   '-G',input$gtffile[[4]],
                   '-p',input$thread,
                   'RNA_seq',input$fastqfile[[4]])
  })
observeEvent(input$runapp,{commandstr <- paste('python3 PCirc_shiny.py -g',input$genomefile[[4]],
                                               '-G',input$gtffile[[4]],
                                               '-p',input$thread,
                                               'RNA_seq',input$fastqfile[[4]])
system(commandstr)})
}



# Run the app ----
shinyApp(ui = ui, server = server)

