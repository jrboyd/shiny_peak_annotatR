#all of the non-reactive server elements
library(shiny)
library(shinyjs)
library(shinyFiles)
library(GenomicRanges)
library(ggplot2)
library(data.table)
library(magrittr)
source("functions_intersect.R")
source("source_gg_venneuler.R")
source("jrb_gg_vennDiagram.R")

bed_path = "~/ShinyApps/shiny_peak_data/beds"
dir.create(bed_path, showWarnings = F)
options(shiny.maxRequestSize=50*1024^2)

# Return the UI for a modal dialog with data selection input. If 'failed' is
# TRUE, then display a message that the previous value was invalid.
dataModal <- function(sets, failed = FALSE) {
  sets_html = HTML(paste(sapply(sets$selected[1], function(x){
    as.character(textInput("TxtRename", label = paste("rename", x), value = x))
  }), collapse = "\n"))
  modalDialog(
    sets_html,
    span('(Please rename the selected sample)'),
    if (failed)
      div(tags$b("One or more names no longer unique!", style = "color: red;")),
    
    footer = tagList(
      modalButton("Cancel"),
      actionButton("BtnConfirmRename", "Confirm")
    )
  )
}

filterModal <- function(failed = FALSE) {
  #unneccessarily complicated for single selection mode but not touching cause it works
  modalDialog(
    span('(Please filter the selected sample)'),
    DT::dataTableOutput("DTPeaksFilter", width = "auto"),
    if (failed)
      div(tags$b("Hey man, this didn't work! Not sure why.", style = "color: red;")),
    
    footer = tagList(
      modalButton("Cancel"),
      actionButton("BtnConfirmFilter", "Confirm")
    ),
    size = "l",
    title = "Filtering"
  )
}



shinyFiles2load = function(shinyF, roots){
  root_path = roots[shinyF$root]
  rel_path = paste0(unlist(shinyF$files), collapse = "/")
  file_path = paste0(root_path, "/", rel_path)
  return(file_path)
}

shinyFiles2save = function(shinyF, roots){
  root_path = roots[shinyF$root]
  rel_path = paste0(unlist(shinyF$name), collapse = "/")
  file_path = paste0(root_path, "/", rel_path)
  return(file_path)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
load_peak_wValidation = function(peak_file, with_notes = F){
  df = read.table(peak_file, stringsAsFactors = F)
  if(ncol(df) == length(peak_cn)){
    if(with_notes){
      showNotification("assuming file is narrowPeak.", type = "warning")
    }else{
      print("assuming file is narrowPeak.")
    }
    colnames(df) = peak_cn  
  }else{
    if(with_notes){
      showNotification("file not narrowPeak, loading as minimal bed file.", type = "warning")
    }else{
      print("file not narrowPeak, loading as minimal bed file.")
    }
    bed_cn = peak_cn[1:4]
    nc = min(ncol(df), length(bed_cn))
    colnames(df)[1:nc] = peak_cn[1:nc]
  }
  return(df)
}

#assumed colnames for MACS2 peak files
peak_cn = c("seqnames", "start", "end", "id", "score", "strand", "FE", "p-value", "q-value", "summit_pos")
#setup root paths
user_roots = dir("/slipstream/home/", full.names = T) %>% paste0(. , "/ShinyData")
user_roots = subset(user_roots, dir.exists(user_roots))
names(user_roots) = dirname(user_roots) %>% basename()
qcframework_load <<- dir("/slipstream/galaxy/uploads/working/qc_framework", pattern = "^output", full.names = T)
names(qcframework_load) <- basename(qcframework_load)
roots_load = c(user_roots, qcframework_load)
roots_output =  c("intersectR" = bed_path, user_roots)