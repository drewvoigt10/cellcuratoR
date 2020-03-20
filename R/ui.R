#' Shiny app server
#'
#' @export
#'
#' @import shiny fluidPage fluidRow column uiOutput tabPanel textOutput h1 h2 h3 h4 h5 br p a imageOutput
#' @importFrom plotly plotlyOutput
#' @importFrom shinydashboard dashboardSidebar sidebarMenu menuItem dashboardBody tabItems tabItem tabBox box dashboardPage dashboardHeader
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinyFiles shinyDirButton
#' @importFrom DT dataTableOutput

sidebar <- dashboardSidebar(
  sidebarMenu(
    id = "menutab",
    menuItem(text = "Dashboard",
             tabName = "dashboard",
             icon = icon("dashboard")
    ),
    shinyDirButton("directory", "Select Seurat object directory", "Please select a folder"),
    uiOutput("dataSelect"),
    menuItem(text = "cellcuratoR",
             tabName = "cellcuratoR",
             icon = icon("file-code-o")
    )
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "dashboard",
            fluidRow(
              tabBox(
                id = "tabset1",
                height = "700px",
                tabPanel("Dimensionality Reduction",
                         textOutput("instructUser"),
                         textOutput("directoryWarning"),
                         plotlyOutput("umapPlot", height = "650px") %>% withSpinner(type = 6, color = "#1F6BFF")
                ),
                tabPanel("Heatmap",
                         plotlyOutput("featurePlot", height = "650px")
                ),
                tabPanel("Differential Expression (Group 1)",
                         plotlyOutput("dePlot1", height = "650px")
                )
              ),
              tabBox(
                id = "tabset2",
                height = "700px",
                tabPanel("Violin Plot",
                         plotOutput("violin_plot", height = "650px")
                ),
                tabPanel("Recluster",
                         plotlyOutput("recluster_plot", height = "650px")
                ),
                tabPanel("Differential Expression (Group 2)",
                         plotlyOutput("dePlot2", height = "650px")
                )
              )
            ),
            fluidRow(
              box(
                column(6,
                       uiOutput("umapHelper"),
                       uiOutput("featurePlotHelper"),
                       uiOutput("deHelper"),
                       uiOutput("deHelper2"),
                       uiOutput("deHelper3"),
                       uiOutput("deHelper5")
                )
              ),
              box(
                column(6,
                       uiOutput("violinHelper"),
                       uiOutput("violinDownload"),
                       uiOutput("reclusterHelper"),
                       uiOutput("reclusterHelper3"),
                       uiOutput("reclusterHelper2"),
                       uiOutput("reclusterHelper4"),
                       uiOutput("reclusterHelper5"),
                       uiOutput("DE_group_1"),
                       uiOutput("deHelper4")
                ),
                column(6,
                       uiOutput("DE_group_2")
                )
              )
            ),
            fluidRow(
              box(
                plotlyOutput("delta_plotly")
              ),
              box(
                plotlyOutput("cluster_visualization")
              )
            ),

            fluidRow(
              DT::dataTableOutput("DE_table")
            )
    ),

    tabItem(tabName = "cellcuratoR",
            h3(strong("Powered by cellcuratoR")),
            fluidRow(
              column(6,
                     br(),
                     br(),
                     p("cellcuratoR is an R package for sharing interactive single-cell expression data from Seurat. Any single-cell RNA
                     sequencing dataset processed with Seurat (v3) can be converted into objects interpretable by cellcuratoR. Code and
                       documentation are available on GitHub, as well as animated instructions for navigating the user interface.", style = "font-size:20px")
              ),
              column(6,
                     imageOutput("cellcuratoR", width = "100%")
              )

            )
    )
  )
)

# Put them together into a dashboardPage
shinyAppUI <- dashboardPage(
  dashboardHeader(title = "cellcuratoR"
  ),
  sidebar,
  body
)
