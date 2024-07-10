##################### SERVER #####################

################################################################################
#Libraries
################################################################################


#Loading Packages
library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(superheat)
library(ggplot2)
#library(gganatogram)
library(BiocManager)
library(plotly)
library(shinycssloaders)
library(shinythemes)
library(RColorBrewer)
library(shinyWidgets)
library(DT)
library(gridExtra)
library(shinyalert)

#options('repos')
#getOption("repos")
# setRepositories()
# options(repos = BiocManager::repositories())

################################################################################
#UIs
################################################################################

First_page<-tabPanel(
  title="Database",
  titlePanel(div(h3("Multi-omics landscape of cancer cell lines", style="margin: 0;"), 
                 h4('v. 2.0 (updated on 25-06-2024)', style="margin: 0;"))),
  #withSpinner(textOutput("lands_1"),color="#0dc5c1"),
  br(),
  withSpinner(plotOutput("lands_1_plot"),color="#0dc5c1"),
  br(),
  titlePanel("Source analyte coverage per Primary Site"),
  br(),
  plotOutput("lands_1_plotly")
)

Data <- tabPanel(
  title = "Analyte Expression",
  titlePanel("Analyte Expression in Cancer Cell lines"),
  sidebarLayout(
    sidebarPanel(
      title = "Expression by Cell-line",
      selectInput("analyte","Choose your module of interest",
                  c("","Gene Expression", "miRNA Expression", "Protein expression",
                    "Metabolite Expression","Drug IC50")),
      selectInput("region","Choose your region of interest",
      c()),
      selectInput("cells","Choose a specific group of cell lines",
                  c("All")),
      selectInput("Histology","Choose a specific histologic group",
                  c()),
      selectInput("Pathology","Choose a specific pathologic group",
                  c()),
      pickerInput("Cells_f","Filter cell lines manually",c(), multiple=T),
      selectInput("Seccion","Choose the top or bottom group",
                  c()),
      actionButton("run_button","Search!",icon=icon("play")),width=3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(
          title = "Analyte profile",
          withSpinner(textOutput("gene_res_plot"),color="#0dc5c1"),
          plotOutput("GeneProf"),
          br(),
          br(),
          br(),
          br(),
          uiOutput("out_resul1")
        ),
        tabPanel(
          title = "Analyte Data",
          textOutput("gene_res"),
          br(),
          br(),
          uiOutput("analyte_resul1"),
          br(),
          dataTableOutput("GeneTable")
        ),
        tabPanel(
          title = "Anatomic expression",
          textOutput("anatomic_exp"),
          textOutput("anatomic_exp2"),
          plotOutput("anatomic_prof"),
          br(),
          br(),
          br(),
          br(),
          br(),
          dataTableOutput("anatomicTable"),
          br(),
          br()
        )
      )
    )
  )
)

ssGSEA <- tabPanel(
  title = "ssGSEA",
  titlePanel("single-sample Gene Set Enrichment Analysis (ssGSEA)"),
  sidebarLayout(
    sidebarPanel(
      title = "Expression by Cell-line",
      selectInput("cells_ss","Filter a specific group of cell lines",
                  c("All")),
      selectInput("Histology_ss","Filter a specific histologic group",
                  c()),
      selectInput("Pathology_ss","Filter a specific pathologic group",
                  c()),
      selectInput("Cells_ss","Choose your cell line manually",c()),
      actionButton("run_button_ss","Search!",icon=icon("play")),width=3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(
          title = "Analyte profile",
          withSpinner(textOutput("gene_res_plot_ss"),color="#0dc5c1"),
          plotOutput("GeneProf_ss"),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          br(),
          uiOutput("out_resul1_ss")
        )
        # ,
        # tabPanel(
        #   title = "Analyte Data",
        #   textOutput("gene_res_ss"),
        #   br(),
        #   br(),
        #   uiOutput("analyte_resul1_ss"),
        #   br(),
        #   dataTableOutput("GeneTable_ss")
        # )
      )
    )
  )
)


Multiple <- tabPanel(
  title = "Multiple Analytes",
  titlePanel("Multiple-Analyte Expression in Cancer Cell lines"),
  sidebarLayout(
    sidebarPanel(
      title = "Expression by Cell-line",
      selectInput("analyte_m","Choose your module of interest",
                  c("","Gene Expression", "miRNA Expression", "Protein expression",
                    "Metabolite Expression")),
      uiOutput("Multiple_cat"),
      textAreaInput("region_m","Write your analytes of interest (one analyte per line)"),
      selectInput("cells_m","Choose a specific group of cell lines",
                  c("All")),
      selectInput("Histology_m","Choose a specific histologic group",
                  c()),
      selectInput("Pathology_m","Choose a specific pathologic group",
                  c()),
      pickerInput("Cells_fm","Filter cell lines manually",c(), multiple=T),
      actionButton("run_button_m","Search!",icon=icon("play")),width=3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(
          title = "Heatmap",
          withSpinner(textOutput("gene_res_plot_m"),color="#0dc5c1"),
          plotOutput("GeneProf_m"),
          br(),
          br(),
          uiOutput("out_resul1_m")
        )
        # ,
        # tabPanel(
        #   title = "Analyte Data",
        #   textOutput("gene_res"),
        #   br(),
        #   br(),
        #   dataTableOutput("GeneTable")
        #)
      )
    )
  )
)



miRgene <- tabPanel(
  title = "miRNA-gene",
  titlePanel("miRNA-gene interactions (based on miRTarBase)"),
  sidebarLayout(
    sidebarPanel(
      title = "miRNA-gene",
      selectInput("mirgene","Filter for (gene or miRNA)",
                  c("","miRNA","gene")),
      selectInput("region1"," ",
                  c()),
      selectInput("region2"," ",
                  c()),
      selectInput("cells2","Choose a specific group of cell lines",
                  c("All")),
      selectInput("Histology2","Choose a specific histologic group",
                  c()),
      selectInput("Pathology2","Choose a specific pathologic group",
                  c()),
      actionButton("run_button2","Search!",icon=icon("play")),width=3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(
          title = "miRNA-gene profile",
          withSpinner(textOutput("miR_gene_text"),color="#0dc5c1"),
          plotlyOutput("miR_gene_plot")
        ),
        tabPanel(
          title = "miRNA-gene Data",
          textOutput("miR_gene_res"),
          br(),
          br(),
          dataTableOutput("miR_geneTable")
        )
      )
    )
  )
)

multi_OmicsR <- tabPanel(
  title = "Multi-omics Data",
  titlePanel("Multi-omics Data"),
  sidebarLayout(
    sidebarPanel(
      title = "Target gene",
      selectizeInput("cell_line_MO","Choose the cell line",
                     c()),
      actionButton("run_button_mO","Search!",icon=icon("play")),width=3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(
          title = "MultiResult",
          dataTableOutput("multiResult_R")
        )
      )
    )
  )
)


silencingR <- tabPanel(
  title = "Modulation tools",
  titlePanel("RNA based tools"),
  sidebarLayout(
    sidebarPanel(
      title = "Target gene",
      selectizeInput("gene_S","Choose the gene",
                  c()),
      actionButton("run_button_s","Search!",icon=icon("play")),width=3
    ),
    
    mainPanel(
       tabsetPanel(
         tabPanel(
            title = "Silencing tools",
           dataTableOutput("silen_DT")
         )
       )
    )
  )
)



ui <- navbarPage(theme = shinytheme("flatly"),
  title = "ShinyThor app", 
  First_page,
  Data,ssGSEA,
  Multiple,
  miRgene,
  #multi_OmicsR,
  silencingR
)




################################################################################
#Required functions
################################################################################

#setwd("/Users/alexis.carrasco/Documents/Shiny ACHILES/")

################################################################################
#Required data
################################################################################
read.csv("Cell_lines_annotations_20181226.txt",sep="\t",fileEncoding="UTF-8")->sample
sample<-sample%>%mutate(name_filt=toupper(gsub(" ","",
                                               gsub("-","",
                                                    gsub("_","",
                                                         gsub("/","",
                                                              gsub("[.]","",
                                                                   gsub("[;]","",
                                                                        gsub("[:]","",
                                                                             gsub("[/]","",
                                                                                  gsub("[)]","",
                                                                                       gsub("[(]","",
                                                                                            gsub("[,]","",Name)))))))))))))


#miRNA
fread("CCLE_miRNA_20180525.gct",encoding = "UTF-8")->miRNA ###Relative Units
miRNA[1:654,]->miRNA
miRNA[,-1]->miRNA
miRNA$Description->names_miRNA
t(miRNA)->miRNA
miRNA[1,]->colnames(miRNA)
miRNA[-1,]->miRNA
miRNA<-data.frame(data.frame(CCLE_ID=rownames(miRNA)),miRNA)
merge(sample[,1:2],miRNA,by="CCLE_ID")->miRNA
miRNA[,-1]->miRNA
colnames(miRNA)<-c("V1",colnames(miRNA)[-1])
gsub('\\.',"-",colnames(miRNA))->colnames(miRNA)

#RNA (gene)
fread("OmicsExpressionProteinCodingGenesTPMLogp1.csv",encoding = "UTF-8")->RNA   ###log2(TPM+1)
sapply(strsplit(colnames(RNA), " "), `[`, 1)->colnames(RNA)

#RPPA Proteins
fread("CCLE_RPPA_20181003.csv",encoding = "UTF-8")->RPPA   ###RPPA Array
#sapply(strsplit(RPPA$V1, "_"), `[`, 1)->RPPA$V1
RPPA <- as.data.frame(RPPA)
colnames(RPPA)[colnames(RPPA) == "V1"] <- "CCLE_ID"
RPPA[,1]->rownames(RPPA)
merge(sample[,1:2],RPPA,by="CCLE_ID")->RPPA
RPPA[,-1]->RPPA
colnames(RPPA)<-c("V1",colnames(RPPA)[-1])

#Metabolites
fread("CCLE_metabolomics_20190502.csv",encoding = "UTF-8")->MTB
MTB <- as.data.frame(MTB)
MTB[,-1]->MTB
colnames(MTB)<-c("V1",colnames(MTB)[-1])

#miRNA-gene
# read_excel(enc2utf8("Data mirna gene.xlsx"))->miRdados
# miRdados<-miRdados%>%mutate(miRNA=gsub("-",".",miRNA))
# unique(miRdados$miRNA)[which(unique(miRdados$miRNA) %in% colnames(miRNA))]->lista_miR
# miRdados[which(miRdados$miRNA %in% lista_miR),]->miRdados
# unique(miRdados$`Target Gene`)[which(unique(miRdados$`Target Gene`) %in% colnames(RNA))]->lista_RNA
# miRdados[which(miRdados$`Target Gene` %in% lista_RNA),]->miRdados
# write.csv(miRdados,"Data_miRNA_gene.csv",fileEncoding = "UTF-8")


#Dados miRTarBase filtrados para incluir regiones humanas y validadas
read.csv("Data_miRNA_gene.csv",fileEncoding = "UTF-8")->miRdados
miRdados[,-1]->miRdados
gsub('\\.',"-",miRdados$miRNA)->miRdados$miRNA

#Data from IC (https://www.cancerrxgene.org/)

read.csv("GDSC_27Oct23.csv",sep=";")->data_IC50
read.csv("GDSC_27Oct23_dic.csv",sep=";")->data_IC50_dic

colnames(data_IC50_dic)<-c("depMapID","COSMIC_ID")


data_IC50%>%
  left_join(data_IC50_dic,"COSMIC_ID")%>%
  mutate(IC_50=exp(LN_IC50))%>%
  left_join(sample%>%select(CCLE_ID,depMapID),"depMapID")%>%
  select(depMapID,DRUG_NAME,IC_50)%>%
  group_by(depMapID,DRUG_NAME)%>%
  summarise(IC_50=max(IC_50))->data_IC50

colnames(data_IC50)<-c("V1","DRUG_NAME","IC_50")
data_IC50%>%
  pivot_wider(id_cols = "V1",names_from="DRUG_NAME",values_from="IC_50")->data_IC50


################################################################################
#Server
################################################################################


server <- function(input, output,session){
  datos <- reactiveValues(x = NULL,ref=NULL,order_color=NULL,
                          final=NULL,o_var=NULL)
  
  # data.frame(`Gene Expression`=sort(colnames(RNA)[-1]),
  #            `miRNA Expression`=c(sort(colnames(miRNA)[-1]),rep(NA,length(colnames(RNA))-length(colnames(miRNA)))),
  #            `Protein Expression`=c(sort(colnames(RPPA)[-1]),rep(NA,length(colnames(RNA))-length(colnames(RPPA)))),
  #            `Metabolite Expression`=c(sort(colnames(MTB)[-1]),rep(NA,length(colnames(RNA))-length(colnames(MTB)))))->availa_m
  # renderDataTable({availa_m
  # })->output$availa_m
  
  
  data.frame(source=c("Gene Expression","miRNA Expression","Protein Expression","Metabolite Expression","Drug IC50"),
             size=c(length(colnames(RNA)[-1]),length(colnames(miRNA)[-1]),
                    length(colnames(RPPA)[-1]),length(colnames(MTB)[-1]),
                    length(colnames(data_IC50)[-1])))->mat_source
  
  data.frame(source=c(rep("Gene expression",nrow(RNA[,1])),
                      rep("Protein expression",length(RPPA[,1])),
                      rep("miRNA expression",length(miRNA[,1])),
                      rep("Metabolite expression",length(MTB[,1])),
                      rep("Drug IC50",nrow(data_IC50[,1]))),
             depMapID=c(as.matrix(RNA[,1]),
                    as.matrix(RPPA[,1]),
                    as.matrix(miRNA[,1]),
                    as.matrix(MTB[,1]),
                    as.matrix(data_IC50[,1])),
             var=1)->tile_graf

  # tile_graf%>%
  #   mutate(var=1)%>%
  #   pivot_wider(names_from=source,values_from=var)%>%
  #   left_join(sample,"depMapID")%>%
  #   mutate(Name=ifelse(is.na(Name),depMapID,paste0(Name," (",depMapID,")")))->tile_graf
    
  
  tile_graf%>%
    left_join(sample,"depMapID")%>%
    mutate(Name=ifelse(is.na(Name),depMapID,paste0(Name," (",depMapID,")")),
           source=factor(source,levels=c("Gene expression","miRNA expression",
                                         "Protein expression","Metabolite expression",
                                         "Drug IC50")))%>%
    group_by(depMapID)%>%
    mutate(n=n())%>%
    arrange(n)->tile_graf
  
  tile_graf%>%
    group_by(Site_Primary)%>%
    mutate(x=n_distinct(depMapID))%>%
    group_by(source,Site_Primary)%>%
    summarise(n=n(),x=unique(x))%>%
    mutate(perc=n/x*100,Site_Primary=ifelse(is.na(Site_Primary),"Others",gsub("[_]"," ",Site_Primary)),
           perc_text=paste0(round(perc,1),"%"))%>%
    ggplot(aes(y=source,x=Site_Primary,size=perc))+geom_point()+theme_bw()+
    geom_text(aes(label=perc_text),color="#996600",size=2.2,vjust=2.8,fontface="bold")+
    theme(axis.text.x=element_text(angle = 45,hjust=1))+ scale_size(range = c(0.00001, 8))+
    scale_size_continuous(name="Percentage \nof cell lines\n covered",
                          breaks=c(25,50,75,100),
                          labels=c("25%","50%","75%","100%"))+xlab("Primary Site")->tile_plotly
  renderPlot({tile_plotly},width=1650,height = 600,res=150)->output$lands_1_plotly

  
  ggplot(mat_source,aes(size,source,fill=source))+geom_bar(stat="identity",alpha=.7)+scale_x_log10()+
    theme_bw()+
    geom_text(aes(label=size), vjust=.4,hjust=1.2,fontface="bold")+
    theme(legend.position = "none"
    )+xlab("Number of analytes")+ scale_fill_brewer(palette = "Set1")->px
  
  
  renderPlot({px},width=800,height = 400,res=150)->output$lands_1_plot
  
  observe({
    
  updateSelectizeInput(session,"gene_S","Choose the gene",
                    c(sort(colnames(RNA)[-1])))
    updateSelectizeInput(session,"cell_line_MO","Choose the cell line",
                         c(sort(sample$name_filt)))
    
    updateSelectInput(session,"cells_ss","Choose a specific group of cell lines",
                      c("All",sort(names(table(sample$Site_Primary)))))
    })
  
  
  ###DATA
  
  observe({
    input$analyte->analyte
    if(analyte==""){
      renderText("")->output$gene_res_plot
      renderPlot({NULL})->output$GeneProf
      output$out_resul1 <- renderUI({NULL})
      updateSelectInput(session,"region","Choose your region of interest",
                        c())
      updateSelectInput(session,"cells","Choose a specific group of cell lines",
                        c())
      updatePickerInput(session,"Cells_f","Filter cell lines manually",c())
      
      updateSelectInput(session,"Seccion","Choose the top or bottom group",
                        c("Top","Bottom"))
      NULL->datos$x
      NULL->datos$ref
    }else if(analyte=="Gene Expression"){
      updateSelectInput(session,"region","Choose your region of interest",
                        c(sort(colnames(RNA[,-1]))))
      updateSelectInput(session,"cells","Choose a specific group of cell lines",
                        c("All",sort(names(table(sample$Site_Primary)))))
      
      updateSelectInput(session,"Seccion","Choose the top or bottom group",
                        c("Top","Bottom"))
      RNA->datos$x
      " log2(TPM+1)"->datos$ref
      renderText("Please, choose your gene")->output$gene_res_plot
      renderPlot({NULL})->output$GeneProf
      output$out_resul1 <- renderUI({NULL})
      
    }else if(analyte=="miRNA Expression"){
      updateSelectInput(session,"region","Choose your region of interest",
                        c(sort(colnames(miRNA[,-1]))))
      updateSelectInput(session,"cells","Choose a specific group of cell lines",
                        c("All",sort(names(table(sample$Site_Primary)))))
      
      updateSelectInput(session,"Seccion","Choose the top or bottom group",
                        c("Top","Bottom"))
      miRNA->datos$x
      " (rel units)"->datos$ref
      renderText("Please, choose your miRNA")->output$gene_res_plot
      renderPlot({NULL})->output$GeneProf
      output$out_resul1 <- renderUI({NULL})
    }else if(analyte=="Protein expression"){
      updateSelectInput(session,"region","Choose your region of interest",
                        c(sort(colnames(RPPA[,-1]))))
      updateSelectInput(session,"cells","Choose a specific group of cell lines",
                        c("All",sort(names(table(sample$Site_Primary)))))
      
      updateSelectInput(session,"Seccion","Choose the top or bottom group",
                        c("Top","Bottom"))
      RPPA->datos$x
      " (rel units)"->datos$ref
      renderText("Please, choose your protein")->output$gene_res_plot
      renderPlot({NULL})->output$GeneProf
      output$out_resul1 <- renderUI({NULL})
      
    }else if(analyte=="Metabolite Expression"){
      updateSelectInput(session,"region","Choose your region of interest",
                        c(sort(colnames(MTB[,-1]))))
      updateSelectInput(session,"cells","Choose a specific group of cell lines",
                        c("All",sort(names(table(sample$Site_Primary)))))
      
      updateSelectInput(session,"Seccion","Choose the top or bottom group",
                        c("Top","Bottom"))
      MTB->datos$x
      " (rel units)"->datos$ref
      renderText("Please, choose your metabolite")->output$gene_res_plot
      renderPlot({NULL})->output$GeneProf
      output$out_resul1 <- renderUI({NULL})
    }else if(analyte=="Drug IC50"){
      updateSelectInput(session,"region","Choose your drug of interest",
                        c(sort(colnames(data_IC50[,-1]))))
      updateSelectInput(session,"cells","Choose a specific group of cell lines",
                        c("All",sort(names(table(sample$Site_Primary)))))
      
      updateSelectInput(session,"Seccion","Choose the top or bottom group",
                        c("Top","Bottom"))
      data_IC50->datos$x
      " (micromoles)"->datos$ref
      renderText("Please, choose your drug")->output$gene_res_plot
      renderPlot({NULL})->output$GeneProf
      output$out_resul1 <- renderUI({NULL})
    }
    
    
    
    # 
    # else if(analyte_m=="Drug IC50"){
    #   
    #   updateSelectInput(session,"cells_m","Choose a specific group of cell lines",
    #                     c("All",sort(names(table(sample$Site_Primary)))))
    #   colnames(data_IC50_dic)<-c("depMapID","COSMIC_ID")
    #   
    #   data_IC50%>%
    #     left_join(data_IC50_dic,"COSMIC_ID")%>%
    #     mutate(IC_50=exp(LN_IC50))->datos$x
    #   " (micromoles)"->datos$ref
    #   renderText("Please, choose your drug")->output$gene_res_plot_m
    #   renderPlot({NULL})->output$GeneProf_m
    #   output$out_resul1_m <- renderUI({NULL})
    # }
    # 
    
  })
  
  observe({
    
    input$cells->cells
    if(cells!="All"){
      
      
      updateSelectInput(session,"Histology",
                        "Choose a specific histologic group",
                        c("All",
                          sort(as.character(sample%>%filter(Site_Primary == cells)%>%select(Histology)%>%distinct()%>%unlist()))))

                        
      updateSelectInput(session,"Pathology","Choose a specific pathologic group",
                        c())
      

  
    } else {
      updateSelectInput(session,"Histology",
                        "Choose a specific histologic group",
                        c("All",
                          sort(as.character(sample%>%select(Histology)%>%distinct()%>%unlist()))))
      
      
      updateSelectInput(session,"Pathology","Choose a specific pathologic group",
                        c())
      
    }
    
    
    
  })
  
  observe({
    input$cells->cells
    input$Histology->Hist
    if(Hist!=""){
    if(Hist!="All"){
      
      
      updateSelectInput(session,"Pathology",
                        "Choose a specific pathologic group",
                        c("All",
                          sort(as.character(sample%>%filter(Site_Primary == cells,
                                                            Histology == Hist)%>%
                                              select(Pathology)%>%distinct()%>%unlist()))))
      

      
    } else {
      
      updateSelectInput(session,"Pathology",
                        "Choose a specific pathologic group",
                        c("All",
                          sort(as.character(sample%>%
                                              select(Pathology)%>%distinct()%>%unlist()))))
      
      
    }
      }
    
    
  })
  
  
  observe({
    input$cells->cells
    input$Histology->Hist
    input$Pathology->Pat
    
    if(cells!="All"){
      sample2<- sample%>%filter(Site_Primary==cells)
    }else{
      sample->sample2
    }
    
    if(Hist!="All"){
      sample2<- sample2%>%filter(Histology==Hist)
    }else{
      sample2->sample2
    }
    
    if(Pat!="All"){
      sample2<- sample2%>%filter(Pathology==Pat)
    }else{
      sample2->sample2
    }
    
    buscar<-sample2%>%
      mutate(Name2=paste0(ifelse(is.na(name_filt),CCLE_ID,name_filt)," (",depMapID,")"))%>%
      select(Name2)%>%unlist()%>%as.character()
    
    updatePickerInput(session,inputId = "Cells_f",label="Filter cell lines manually",
                      choices=buscar,
                      selected=buscar,
                      #multiple = TRUE,
                      options = list(`actions-box` = TRUE))
  })
  
  observeEvent(input$run_button,{
    
    input$cells->cells
    input$Histology->Hist
    input$Pathology->Pat
    input$region->gene
    input$Seccion->seccion
    datos$x->specificSource
    input$Cells_f->cells_F
    input$analyte->analyte
    
    ##### Create anatogram
    
    # specificSource%>%
    #   select(c("V1",gene))%>%as.data.frame()->resul_anato
    # resul_anato<-merge(sample%>%select(depMapID,Name,Site_Primary,
    #                                Histology,Pathology,Gender),resul_anato,by.y="V1",by.x="depMapID")
    # 
    # resul_anato<-resul_anato%>%filter(Pathology=="primary", Gender %in% c("male","female"))
    # 
    # resul_anato[,c(3,6,7)]->vale_anato
    # colnames(vale_anato)<-c("organ","Gender","no_scale_value")
    # 
    # 
    # vale_anato<-vale_anato%>%
    #   mutate(no_scale_value=as.numeric(no_scale_value))
    # 
    # gsub("central_nervous_system","frontal_cortex",vale_anato$organ)->vale_anato$organ
    # gsub("autonomic_ganglia","spinal_cord",vale_anato$organ)->vale_anato$organ
    # gsub("biliary_tract","gall_bladder",vale_anato$organ)->vale_anato$organ
    # gsub("haematopoietic_and_lymphoid_tissue","lymph_node",vale_anato$organ)->vale_anato$organ
    # gsub("large_intestine","colon",vale_anato$organ)->vale_anato$organ
    # gsub("oesophagus","esophagus",vale_anato$organ)->vale_anato$organ
    # gsub("soft_tissue","nerve",vale_anato$organ)->vale_anato$organ
    # gsub("thyroid","thyroid_gland",vale_anato$organ)->vale_anato$organ
    # gsub("upper_aerodigestive_tract","trachea",vale_anato$organ)->vale_anato$organ
    # 
    # vale_anato<-vale_anato%>%
    #       dplyr::filter( organ %in% unique(hgMale_key$organ,hgFemale_key$organ))%>%
    #       mutate(value=scale(as.numeric(no_scale_value)))
    # 
    # 
    # vale_anato%>%
    #   filter(Gender=="male")%>%
    #   inner_join(hgMale_key[,-4],by = "organ")%>%
    #   group_by(organ,colour)%>%
    #   reframe(value=mean(na.omit(as.numeric(value))))%>%
    #   arrange(value)->vale_anato_h
    # 
    # 
    # vale_anato%>%
    #   filter(Gender=="female")%>%
    #   inner_join(hgFemale_key[,-4],by = "organ")%>%
    #   group_by(organ,colour)%>%
    #   reframe(value=mean(na.omit(as.numeric(value))))%>%
    #   arrange(value)->vale_anato_f
    # 
    # 
    # anato_h<-gganatogram(data=vale_anato_h, fillOutline='white', organism='human', sex='male', fill="value")+ 
    #   theme_void() +
    #   scale_fill_gradient(low = "blue", high = "red",name="Scaled \nvalues") 
    # 
    # anato_f<-gganatogram(data=vale_anato_f, fillOutline='white', organism='human', sex='female', fill="value")+ 
    #   theme_void() +
    #   scale_fill_gradient(low = "blue", high = "red",name="Scaled \nvalues") 
    # 
    # 
    # renderText({paste0("Anatomic levels of ",gene," ",datos$ref," . Mean expression per organ.")})->output$anatomic_exp
    # renderText({"Cells derived from male or female samples are shown in left or right panels, respectively"})->output$anatomic_exp2
    # renderPlot({
    #   grid.arrange(anato_h, anato_f, ncol=2)->anatomic_prof
    #   
    #   print(anatomic_prof)},width=1000,height = 500,res=150)->output$anatomic_prof
    # 
    # vale_anato<-vale_anato%>%
    #   group_by(organ,Gender)%>%
    #   reframe(no_scale_value=mean(na.omit(as.numeric(no_scale_value))))
    # 
    # vale_anato$no_scale_value<-round(as.numeric(vale_anato$no_scale_value),2)
    # vale_anato<-vale_anato%>%arrange(-no_scale_value)
    # colnames(vale_anato)<-c("Affected Organ","Gender",paste0("Expression - ",datos$ref))
    # 
    # renderDataTable({vale_anato})->output$anatomicTable
    # 
    #### End anatogram
    
    if(cells!="All"){
      sample2<- sample%>%filter(Site_Primary==cells)
    }else{
      sample->sample2
    }
    
    if(Hist!="All"){
      sample2<- sample2%>%filter(Histology==Hist)
    }else{
      sample2->sample2
    }
    
    if(Pat!="All"){
      sample2<- sample2%>%filter(Pathology==Pat)
    }else{
      sample2->sample2
    }
    
    sample2%>%
      mutate(Name2=paste0(ifelse(is.na(name_filt),CCLE_ID,name_filt)," (",depMapID,")"))->sample2
    
    buscar<-sample2%>%filter(Name2 %in% cells_F)%>%select(depMapID)%>%unlist()%>%as.character()
    
    result<-specificSource%>%filter(V1 %in% buscar)%>%
      select(c("V1",gene))%>%as.data.frame()
    
    result<-merge(sample2%>%select(depMapID,Name,Site_Primary,
                                  Histology,Pathology),result,by.y="V1",by.x="depMapID")

    colnames(result)<-c("DepMap_ID","Name","Site","Histology","Pathology","gene")
    
    
    if(seccion=="Top"){
      numbers<-10
    }else if (seccion=="Bottom"){
      numbers<- 10*-1
    }
    p1<-result%>%
      mutate(gene=as.numeric(gene))%>%
      arrange(gene)%>%
      top_n(numbers)%>%
      mutate(Name = fct_reorder(Name, desc(gene)))%>%
      ggplot(aes(Name,gene))+geom_col(fill="#1B9E77")+
      ylab(paste0(gene,datos$ref))+
      xlab("Cell line")+
      theme_bw()+ 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    renderText({paste0(seccion, " ",gene,ifelse(analyte=="Drug IC50"," IC50 levels in "," expression in "),
                       cells," cell lines (Histology: ",
                       Hist,", Pathology: ",Pat,")")})->output$gene_res_plot
    
    renderPlot({p1},width = 145+ifelse(nrow(result)>15,530,nrow(result)*37),height = 450, res=150)->output$GeneProf
    
    
    plotdCt<-reactive({p1})
    output$out_resul1 <- renderUI({
      
      downloadButton("downloadPlotdCt","Download this plot")         
    })
    
    output$downloadPlotdCt <- downloadHandler(
      filename = paste0(seccion, " ",gene," expression in ",cells," (H_", Hist,"_P_",Pat, ").png"),
      content = function(file) {
        png(file,width = 540+ifelse(nrow(result)>15,2160,nrow(result)*154),height = 1800, res=600)
        print(plotdCt())
        dev.off()
      }) 
    
    
    
    
    result<-result%>%mutate(gene=round(as.numeric(gene),2))%>%arrange(-gene)
    
    colnames(result)<-c("DepMap ID","Cell line","Primary Site","Histology","Pathology",paste0(gene,datos$ref))
  
    renderText({paste0("After your selection, we filtered ",nrow(result)," out of 1461 cell lines")})->output$gene_res
    
    output$analyte_resul1 <- renderUI({
      
      downloadButton("downloadConvertedData","Download this table")         
    })
    
    
    output$downloadConvertedData <- downloadHandler(
      filename = function() {
        gsub(":","_",Sys.time())->timed
        paste("Analyte_Data_",gene,"_", timed, ".csv", sep="")
      },
      content = function(file) {
        write.csv(result, file,row.names = F)
      }
    )
    
    renderDataTable({result})->output$GeneTable
    
    
    })
  
  
  ### ssGSEA
  

  observe({
    renderText("")->output$gene_res_plot_ss
    
    input$cells_ss->cells_ss
    if(cells_ss!="All"){
      
      
      updateSelectInput(session,"Histology_ss",
                        "Filter a specific histologic group",
                        c("All",
                          sort(as.character(sample%>%filter(Site_Primary == cells_ss)%>%
                                              select(Histology)%>%distinct()%>%unlist()))))
      
      
      updateSelectInput(session,"Pathology_ss","Filter a specific pathologic group",
                        c())
      
      
      
    } else {
      updateSelectInput(session,"Histology_ss",
                        "Filter a specific histologic group",
                        c("All",
                          sort(as.character(sample%>%select(Histology)%>%distinct()%>%unlist()))))
      
      
      updateSelectInput(session,"Pathology_ss","Filter a specific pathologic group",
                        c())
      
    }
    
    
    
  })
  
  observe({
    input$cells_ss->cells_ss
    input$Histology_ss->Hist_ss
    if(Hist_ss!=""){
      if(Hist_ss!="All"){
        
        
        updateSelectInput(session,"Pathology_ss",
                          "Filter a specific pathologic group",
                          c("All",
                            sort(as.character(sample%>%filter(Site_Primary == cells_ss,
                                                              Histology == Hist_ss)%>%
                                                select(Pathology)%>%distinct()%>%unlist()))))
        
        
        
      } else {
        
        updateSelectInput(session,"Pathology_ss",
                          "Filter a specific pathologic group",
                          c("All",
                            sort(as.character(sample%>%
                                                select(Pathology)%>%distinct()%>%unlist()))))
        
        
      }
    }
    
    
  })
  
  
  observe({
    input$cells_ss->cells_ss
    input$Histology_ss->Hist_ss
    input$Pathology_ss->Pat_ss
    
    if(cells_ss!="All"){
      sample2<- sample%>%filter(Site_Primary==cells_ss)
    }else{
      sample->sample2
    }
    
    if(Hist_ss!="All"){
      sample2<- sample2%>%filter(Histology==Hist_ss)
    }else{
      sample2->sample2
    }
    
    if(Pat_ss!="All"){
      sample2<- sample2%>%filter(Pathology==Pat_ss)
    }else{
      sample2->sample2
    }
    
    buscar<-sample2%>%
      mutate(Name2=paste0(ifelse(is.na(name_filt),CCLE_ID,name_filt)," (",depMapID,")"))%>%
      select(Name2)%>%unlist()%>%as.character()
    
    updateSelectInput(session,inputId = "Cells_ss",label="Filter cell lines manually",
                      buscar)
  })
  
  observeEvent(input$run_button_ss,{
    renderText({paste0("Loading...")})->output$gene_res_plot_ss
    
    input$cells_ss->cells_ss
    input$Histology_ss->Hist_ss
    input$Pathology_ss->Pat_ss
    datos$x->specificSource
    input$Cells_ss->Cells_ss
  
    str_remove(Cells_ss, ".*[(]")->Cells_ss
    substr(Cells_ss,1,10)->Cells_ss
    
    sample%>%
      select(depMapID,Name,Pathology,Site_Primary,Histology)->samp2
    
    colnames(samp2)<-c("SAMPLE_ID","Name","Pathology","Site_Primary","Histology")
    
    # ssGSEA GeneOntology
    ssGSEA_CellLines_C7_immune <- fread("ssGSEA_CellLines_C5_gene_ontology.txt")
    
    samp2%>%
      filter(SAMPLE_ID==Cells_ss)%>%
      left_join(ssGSEA_CellLines_C7_immune%>%select(-V1),"SAMPLE_ID")->samp3
    
    
    samp3%>%
      filter(SAMPLE_ID==Cells_ss)%>%
      pivot_longer(!c("SAMPLE_ID","Name","Pathology","Site_Primary","Histology"),
                   names_to="Pathway",values_to="Score")%>%
      mutate(Score=as.numeric(Score),Score_sig=ifelse(Score>0,"Pos","Neg"),
             Pathway=gsub("[_]"," ",Pathway),
             Pathway=sapply(Pathway, function(x) paste(strwrap(x, 40), collapse = "\n")))%>%
      filter(abs(Score)>0.1)%>%
      group_by(Score_sig)%>%
      top_n(10,abs(Score))%>%
      ungroup()%>%
      arrange(Score)%>%
      mutate(Pathway = fct_reorder(Pathway, Score))%>%
      ggplot(aes(x=Score,y=Pathway,color=Score_sig))+
      geom_segment(aes(x=0,xend=Score,size=.2), color="grey40")+geom_point(aes(size=abs(Score)))+
      scale_size_continuous(guide="none")+theme_bw()+ 
      scale_color_brewer(guide="none",palette = "Dark2")->p1
    

     renderText({paste0("Relevant pathways enriched in the ", Cells_ss, " cell line (SSGSEA).")})->output$gene_res_plot_ss
    # 
     renderPlot({p1},width = 800,height = 850, res=120)->output$GeneProf_ss
    # 
    # 
     plot_ssGSEA<-reactive({p1})
    output$out_resul1_ss <- renderUI({

      downloadButton("downloadssGSEA","Download this plot")
    })

    output$downloadssGSEA <- downloadHandler(
      filename = paste0("Relevant pathways enriched in the ", Cells_ss, " cell line (SSGSEA).png"),
      content = function(file) {
        png(file,width = 3200,height = 3400, res=600)
        print(plot_ssGSEA())
        dev.off()
      })
    # 
    # 
    # 
    # 
    # result<-result%>%mutate(gene=round(as.numeric(gene),2))%>%arrange(-gene)
    # 
    # colnames(result)<-c("DepMap ID","Cell line","Primary Site","Histology","Pathology",paste0(gene,datos$ref))
    # 
    # renderText({paste0("After your selection, we filtered ",nrow(result)," out of 1461 cell lines")})->output$gene_res
    # 
    # output$analyte_resul1 <- renderUI({
    #   
    #   downloadButton("downloadConvertedData","Download this table")         
    # })
    # 
    # 
    # output$downloadConvertedData <- downloadHandler(
    #   filename = function() {
    #     gsub(":","_",Sys.time())->timed
    #     paste("Analyte_Data_",gene,"_", timed, ".csv", sep="")
    #   },
    #   content = function(file) {
    #     write.csv(result, file,row.names = F)
    #   }
    # )
    # 
    # renderDataTable({result})->output$GeneTable
    
    
  })
  
  
  
  ###MULTIPLE

  observe({
    input$analyte_m->analyte_m
    
    if(analyte_m!=""){
    output$Multiple_cat<-
      renderUI({
        
        actionButton("Multiple_cat_bttn","Check available analytes")
      })
    }
  })
  
  observe({
    input$analyte_m->analyte_m
    if(analyte_m==""){
      renderText("")->output$gene_res_plot_m
      renderPlot({NULL})->output$GeneProf_m
      output$out_resul1_m <- renderUI({NULL})
      updateSelectInput(session,"cells_m","Choose a specific group of cell lines",
                        c())
      updatePickerInput(session,"Cells_fm","Filter cell lines manually",c())
      
      NULL->datos$x
      NULL->datos$ref
    }else if(analyte_m=="Gene Expression"){

      updateSelectInput(session,"cells_m","Choose a specific group of cell lines",
                        c("All",sort(names(table(sample$Site_Primary)))))

      RNA->datos$x
      " log2(TPM+1)"->datos$ref
      renderText("Please, write your genes")->output$gene_res_plot_m
      renderPlot({NULL})->output$GeneProf_m
      output$out_resul1_m <- renderUI({NULL})
      

      
    }else if(analyte_m=="miRNA Expression"){

      updateSelectInput(session,"cells_m","Choose a specific group of cell lines",
                        c("All",sort(names(table(sample$Site_Primary)))))

      miRNA->datos$x
      " (rel units)"->datos$ref
      renderText("Please, choose your miRNAs")->output$gene_res_plot_m
      renderPlot({NULL})->output$GeneProf_m
      output$out_resul1_m <- renderUI({NULL})
    }else if(analyte_m=="Protein expression"){

      updateSelectInput(session,"cells_m","Choose a specific group of cell lines",
                        c("All",sort(names(table(sample$Site_Primary)))))
      
      RPPA->datos$x
      " (rel units)"->datos$ref
      renderText("Please, choose your proteins")->output$gene_res_plot_m
      renderPlot({NULL})->output$GeneProf_m
      output$out_resul1_m <- renderUI({NULL})
      
    }else if(analyte_m=="Metabolite Expression"){

      updateSelectInput(session,"cells_m","Choose a specific group of cell lines",
                        c("All",sort(names(table(sample$Site_Primary)))))

      MTB->datos$x
      " (rel units)"->datos$ref
      renderText("Please, choose your metabolite")->output$gene_res_plot_m
      renderPlot({NULL})->output$GeneProf_m
      output$out_resul1_m <- renderUI({NULL})
    }else if(analyte_m=="Drug IC50"){
      
      updateSelectInput(session,"cells_m","Choose a specific group of cell lines",
                        c("All",sort(names(table(sample$Site_Primary)))))
      
      data_IC50->datos$x
      " (micromoles)"->datos$ref
      renderText("Please, choose your drug")->output$gene_res_plot_m
      renderPlot({NULL})->output$GeneProf_m
      output$out_resul1_m <- renderUI({NULL})
    }
    
  })
  
  observeEvent(input$Multiple_cat_bttn, {
    datos$x->specificSource
    input$analyte_m->analyte_m
    
    shinyalert(
      html = TRUE,
      text = tagList(
        selectInput("available_analytes", analyte_m, sort(colnames(specificSource[,-1])))
      ),
      closeOnClickOutside=T
    )
  })
  
  
  observe({
    
    input$cells_m->cells_m
    if(cells_m!="All"){
      
      
      updateSelectInput(session,"Histology_m",
                        "Choose a specific histologic group",
                        c("All",
                          sort(as.character(sample%>%filter(Site_Primary == cells_m)%>%select(Histology)%>%distinct()%>%unlist()))))
      
      
      updateSelectInput(session,"Pathology_m","Choose a specific pathologic group",
                        c())
      
    } else {
      updateSelectInput(session,"Histology_m",
                        "Choose a specific histologic group",
                        c("All",
                          sort(as.character(sample%>%select(Histology)%>%distinct()%>%unlist()))))
      
      
      updateSelectInput(session,"Pathology_m","Choose a specific pathologic group",
                        c())
      
    }
    
    
    
  })
  
  observe({
    input$cells_m->cells_m
    input$Histology_m->Hist_m
    if(Hist_m!=""){
      if(Hist_m!="All"){
        
        
        updateSelectInput(session,"Pathology_m",
                          "Choose a specific pathologic group",
                          c("All",
                            sort(as.character(sample%>%filter(Site_Primary == cells_m,
                                                              Histology == Hist_m)%>%
                                                select(Pathology)%>%distinct()%>%unlist()))))
        
        
        
      } else {
        
        updateSelectInput(session,"Pathology_m",
                          "Choose a specific pathologic group",
                          c("All",
                            sort(as.character(sample%>%
                                                select(Pathology)%>%distinct()%>%unlist()))))
        
        
      }
    }
    
    
  })

  observe({
    input$cells_m->cells_m
    input$Histology_m->Hist_m
    input$Pathology_m->Pat_m
    
    if(cells_m!="All"){
      sample2<- sample%>%filter(Site_Primary==cells_m)
    }else{
      sample->sample2
    }
    
    if(Hist_m!="All"){
      sample2<- sample2%>%filter(Histology==Hist_m)
    }else{
      sample2->sample2
    }
    
    if(Pat_m!="All"){
      sample2<- sample2%>%filter(Pathology==Pat_m)
    }else{
      sample2->sample2
    }
    
    buscar<-sample2%>%
      mutate(Name2=paste0(ifelse(is.na(name_filt),CCLE_ID,name_filt)," (",depMapID,")"))%>%
      select(Name2)%>%unlist()%>%as.character()
    
    updatePickerInput(session,inputId = "Cells_fm",label="Filter cell lines manually",
                      choices=buscar,
                      selected=buscar,
                      #multiple = TRUE,
                      options = list(`actions-box` = TRUE))
  })
  
  observeEvent(input$run_button_m,{
    input$cells_m->cells_m
    input$Histology_m->Hist_m
    input$Pathology_m->Pat_m
    input$region_m->gene_m
    datos$x->specificSource
    input$Cells_fm->cells_Fm
    
    if(cells_m!="All"){
      sample2<- sample%>%filter(Site_Primary==cells_m)
    }else{
      sample->sample2
    }
    
    if(Hist_m!="All"){
      sample2<- sample2%>%filter(Histology==Hist_m)
    }else{
      sample2->sample2
    }
    
    if(Pat_m!="All"){
      sample2<- sample2%>%filter(Pathology==Pat_m)
    }else{
      sample2->sample2
    }
    
    sample2%>%
      mutate(Name2=paste0(ifelse(is.na(name_filt),CCLE_ID,name_filt)," (",depMapID,")"))->sample2
    
    buscar<-sample2%>%filter(Name2 %in% cells_Fm)%>%select(depMapID)%>%unlist()%>%as.character()
    

    if(cells_m!="All"){
    str_split(gene_m," ")->gene_m
    str_split(gene_m[[1]],"\n")->gene_m
    unlist(gene_m)->gene_m
    str_trim(gene_m)->gene_m
    gene_m [which (gene_m %in% colnames(specificSource))]->gene_m
  
    if(length(gene_m)>1){
      result<-specificSource%>%filter(V1 %in% buscar)%>%
        select(c("V1",gene_m))%>%as.data.frame()
      
      result<-merge(sample2%>%select(depMapID,Name,Site_Primary,
                                     Histology,Pathology),result,by.y="V1",by.x="depMapID")
      
      colnames(result)<-c("DepMap_ID","Name","Site","Histology","Pathology",colnames(result)[6:ncol(result)])
      
      apply(result[,6:ncol(result)],2,as.numeric)->result[,6:ncol(result)]
      
      result$Name->rownames(result)
      
      colMain <- colorRampPalette(brewer.pal(8, "Blues"))(25)
      
      renderText({paste0("Expression of ",length(gene_m)," regions in ", cells_m," cell lines (Histology: ",
                         Hist_m,", Pathology: ",Pat_m,")")})->output$gene_res_plot_m
      
      
      library(ComplexHeatmap)
      
      
      length(unique(result$Pathology))->Patho
      length(unique(result$Histology))->Histo
      
      if(Patho>1 & Histo >1 ){
        ha <- HeatmapAnnotation(
          Pathology = result$Pathology, Histology = result$Histology)
        
        P_M1<-Heatmap(as.matrix(t(result[,6:ncol(result)])), name = "Expression", 
                      top_annotation = ha,col = colorRampPalette(brewer.pal(8, "Blues"))(25))
      }else if (Patho>1){
        ha <- HeatmapAnnotation(
          Pathology = result$Pathology)
        
        P_M1<-Heatmap(as.matrix(t(result[,6:ncol(result)])), name = "Expression", 
                      top_annotation = ha,col = colorRampPalette(brewer.pal(8, "Blues"))(25))
      }else if(Histo >1){
        ha <- HeatmapAnnotation(
          Histology = result$Histology)
        
        P_M1<-Heatmap(as.matrix(t(result[,6:ncol(result)])), name = "Expression", 
                      top_annotation = ha,col = colorRampPalette(brewer.pal(8, "Blues"))(25))
      }else{
        
        P_M1<-Heatmap(as.matrix(t(result[,6:ncol(result)])), name = "Expression",
                      col = colorRampPalette(brewer.pal(8, "Blues"))(25))
      }

      

      
      
      renderPlot({
        
        P_M1
      },
      width = 120+nrow(result)*40,
      height = 240+length(gene_m)*30, res=100)->output$GeneProf_m
      
      
      
      plotdCt_m<-reactive({P_M1})
      output$out_resul1_m <- renderUI({
        
        downloadButton("downloadPlotdCt_m","Download this plot")         
      })
      
      output$downloadPlotdCt_m <- downloadHandler(
        filename = paste0("Expression of ",length(gene_m),"genes in ", cells_m," cell lines (Histology: ",
                          Hist_m,", Pathology: ",Pat_m,")",".png"),
        content = function(file) {
          png(file,width = 1000+nrow(result)*125,
              height = 1600+length(gene_m)*200, res=600)
          print(plotdCt_m())
          dev.off()
        }) 
    }else{
      renderText("Please, include more than one analyte to be analyzed")->output$gene_res_plot_m
    }
    
    }
    
    
    # 
    # result<-result%>%mutate(gene=as.numeric(gene))%>%arrange(-gene)
    # 
    # colnames(result)<-c("DepMap ID","Cell line","Primary Site","Histology","Pathology",paste0(gene,datos$ref))
    # 
    # renderText({paste0("After your selection, we filtered ",nrow(result)," out of 1461 cell lines")})->output$gene_res
    # 
    # 
    # 
    # renderDataTable({result})->output$GeneTable
    # 
    
  })
  
  
  ### miRgene
  
  observe({
    input$mirgene->mirgene
    if(mirgene=="miRNA"){
      updateSelectInput(session,"region1","Choose your miRNA of interest",
                        c(sort(unique(miRdados$miRNA))))
      updateSelectInput(session,"region2","Choose your gene of interest",
                        c("Choose your miRNA first"))
      updateSelectInput(session,"cells2","Choose a specific group of cell lines",
                        c("All",sort(names(table(sample$Site_Primary)))))
      

    }else if(mirgene=="gene"){
      updateSelectInput(session,"region1","Choose your gene of interest",
                        c(sort(unique(miRdados$Target.Gene))))
      updateSelectInput(session,"region2","Choose your miRNA of interest",
                        c("Choose your gene first"))
      updateSelectInput(session,"cells2","Choose a specific group of cell lines",
                        c("All",sort(names(table(sample$Site_Primary)))))
      
    }else if(mirgene==""){
      updateSelectInput(session,"region1","Choose your region of interest",
                        c())
      updateSelectInput(session,"cells2","Choose a specific group of cell lines",
                        c())

    }
    
  })
  
  observe({
    input$mirgene->mirgene
    input$region1->region1
    
    if(mirgene != ""){
      if(mirgene == "miRNA"){
        updateSelectInput(session,"region2","Choose your gene of interest",
                          c(sort(unique(miRdados[which(miRdados$miRNA == region1),]$Target.Gene))))
        renderText("Please, choose your miRNA")->output$miR_gene_text
      }else if (mirgene == "gene"){
        updateSelectInput(session,"region2","Choose your miRNA of interest",
                          c(sort(unique(miRdados[which(miRdados$Target.Gene == region1),]$miRNA))))
        renderText("Please, choose your gene")->output$miR_gene_text
      }
    }else{
      renderText("")->output$miR_gene_text
    }
    

    
  })
  
  observe({
    
    input$cells2->cells2
    if(cells2!="All"){
      
      
      updateSelectInput(session,"Histology2",
                        "Choose a specific histologic group",
                        c("All",
                          sort(as.character(sample%>%
                                              filter(Site_Primary == cells2)%>%
                                              select(Histology)%>%
                                              distinct()%>%unlist()))))
      
      
      updateSelectInput(session,"Pathology2","Choose a specific pathologic group",
                        c())
      
    } else {
      updateSelectInput(session,"Histology2",
                        "Choose a specific histologic group",
                        c("All",
                          sort(as.character(sample%>%select(Histology)%>%
                                              distinct()%>%unlist()))))
      
      
      updateSelectInput(session,"Pathology2","Choose a specific pathologic group",
                        c())
      
    }
    
    
    
  })
  
  observe({
    input$cells2->cells2
    input$Histology2->Hist2
    if(Hist2!=""){
      if(Hist2!="All"){
        
        
        updateSelectInput(session,"Pathology2",
                          "Choose a specific pathologic group",
                          c("All",
                            sort(as.character(sample%>%filter(Site_Primary == cells2,
                                                              Histology == Hist2)%>%
                                                select(Pathology)%>%distinct()%>%unlist()))))
        
        
        
      } else {
        
        updateSelectInput(session,"Pathology2",
                          "Choose a specific pathologic group",
                          c("All",
                            sort(as.character(sample%>%
                                                select(Pathology)%>%distinct()%>%unlist()))))
        
        
      }
    }
    
    
  })
  
  observeEvent(input$run_button2,{
    input$cells2->cells2
    input$Histology2->Hist2
    input$Pathology2->Pat2
    input$mirgene->mirgene
    input$region1->reg1
    input$region2->reg2
    renderText(" ")->output$miR_gene_text
    
    if(cells2!="All"){
      sample2<- sample%>%filter(Site_Primary==cells2)
    }else{
      sample->sample2
    }
    
    if(Hist2!="All"){
      sample2<- sample2%>%filter(Histology==Hist2)
    }else{
      sample2->sample2
    }
    
    if(Pat2!="All"){
      sample2<- sample2%>%filter(Pathology==Pat2)
    }else{
      sample2->sample2
    }
    
    buscar<-sample2%>%select(depMapID)%>%unlist()%>%as.character()
  
    if(mirgene != ""){
      if(mirgene =="miRNA"){
        tr1<-miRNA%>%select(all_of(c("V1",reg1)))
        tr2<-RNA%>%select(all_of(c("V1",reg2)))
        merge(tr1,tr2,by="V1")->tr
        2^tr[,3]-1-> tr[,3]
        tr<-merge(sample2%>%select(depMapID,Name,Site_Primary,
                                   Histology,Pathology),tr,by.x="depMapID",by.y="V1")
        reg1->miR_name
        reg2->gene_name
      }else if(mirgene =="gene"){
        tr1<-miRNA%>%select(all_of(c("V1",reg2)))
        tr2<-RNA%>%select(all_of(c("V1",reg1)))
        merge(tr1,tr2,by="V1")->tr
        2^tr[,3]-1-> tr[,3]
        tr<-merge(sample2%>%select(depMapID,Name,Site_Primary,
                               Histology,Pathology),tr,by.x="depMapID",by.y="V1")
        reg2->miR_name
        reg1->gene_name
      }
      
      colnames(tr)<-c("DepMap_ID","Name","Site","Histology","Pathology","miRNA","gene")
      
      tr<-tr%>%mutate(miRNA=round(as.numeric(miRNA),2),
                      gene=round(as.numeric(gene),2))
      
      colnames(tr)<-c("DepMap_ID","Name","Site","Histology","Pathology",miR_name,gene_name)
      
      fig<-style(ggplot(tr,aes(get(miR_name),get(gene_name)))+
                      geom_point(color="darkblue",alpha=.7)+
                      geom_smooth(method="lm",alpha=0.1,color="orange")+
                      theme_bw()+
                      xlab(paste0(miR_name," (rel units)"))+
                      ylab(paste0(gene_name," (TPM)")),text=paste0(tr$Name," (",tr$DepMap_ID,")\n",
                                                  miR_name,": ",round(tr[,miR_name],2),"\n",
                                                  gene_name,": ",round(tr[,gene_name],2)))|>
        add_annotations(
          xref = "paper", yref = "paper",
          x = .95 , y = .95, 
          text = paste0("<i>R</i> = ", round(cor(tr[,miR_name], tr[,gene_name]),2),
                        "<br>",
                        "<i>P</i> = ", formatC(cor.test(tr[,miR_name], tr[,gene_name])$p.value,
                                               format="e", digits=2)),
          showarrow = F,  
          align = "left") 
      
      output$miR_gene_plot<-renderPlotly(fig)
      renderDataTable({tr},
                      options = list(
                        autoWidth = TRUE))->output$miR_geneTable
    }
    
    
    
  })

  ### multi-omics
  
  observeEvent(input$run_button_mO,{
    input$cell_line_MO->cell_line_mO
    
    sample[which(sample$name_filt==cell_line_mO),]$depMapID->ID
    
   head(t(as.matrix(RNA[which(RNA$V1==ID),-1])))
  })
  
  ### silencing
  
  observeEvent(input$run_button_s,{
    
    input$gene_S->gene_S
    
    if(gene_S!=""){
      
      data.frame(`Target Gene`=c("miRNA (miRTarBase)","cirRNA (circInteractome)"),
                 gene_S=c(paste0(unique(miRdados[which(miRdados$Target.Gene == gene_S),]$miRNA),collapse = ", "),
                          paste0("<a href='  https://circinteractome.nia.nih.gov/api/v2/circsearch?circular_rna_query=&gene_symbol_query=",
                          gene_S,"&submit=circRNA+Search' target='_blank'>",gene_S,"</a>"))
                 )->silencing
      
      colnames(silencing)<-c("Target Gene",gene_S)
      
      renderDataTable({datatable(silencing, 
                                 options = list(pageLength = 15, lengthChange = FALSE),
                                 escape = FALSE,rownames=F)})->output$silen_DT
      
    }
   
    })
  
  }

shinyApp(ui = ui, server = server)
