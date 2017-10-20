
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyjs)
library(shinyFiles)
library(GenomicRanges)
library(ggplot2)
library(data.table)
library(magrittr)
source("functions_intersect.R")
source("source_gg_venneuler.R")

bed_path = "~/ShinyApps/shiny_peak_data/beds"

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

shinyServer(function(input, output, session) {
  user_roots = dir("/slipstream/home/", full.names = T) %>% dir(. ,pattern = "^ShinyData$", full.names = T)
  names(user_roots) = dirname(user_roots) %>% basename()
  qcframework_load <<- dir("/slipstream/galaxy/uploads/working/qc_framework", pattern = "^output", full.names = T)
  names(qcframework_load) <- basename(qcframework_load)
  roots_load = c(user_roots, qcframework_load)
  shinyFileChoose(input, 'FilesLoadData', roots= roots_load, filetypes=c("narrowPeak", "broadPeak"))
  roots_output =  c("bed_dir" = bed_path)
  dir.create(roots_output, showWarnings = F)
  shinyFileSave(input, 'FilesSaveResults', roots= roots_output, filetypes=c("bed"))
  peak_cn = c("seqnames", "start", "end", "id", "score", "strand", "FE", "p-value", "q-value", "summit_pos")
  peak_prev_file = reactiveVal(value = "", label = "peak_prev_file")
  peak_prev_name = reactiveVal(value = "", label = "peak_prev_name")
  #stores all set data, keyed by set_id
  peak_files = reactiveVal(value = list(), label = "peak_files")
  #stores the order of peak_files as well as any updated names
  
  observeEvent(input$BtnAddFile, {
    if(is.null(peak_dfr())) return(NULL)
    tmp = peak_files()
    tmp[[input$TxtFileName]] = peak_dfr()
    peak_files(tmp)
    sets = input$ChooseIntervalSets
    output$SetChooser = renderUI({
      pfl = sets$left
      #add to active
      pfr = c(sets$right, isolate(input$TxtFileName))
      
      # print(pf)
      return(chooserInput(inputId = "ChooseIntervalSets", 
                          leftLabel = "Ready", rightLabel = "For Analysis", 
                          leftChoices = pfl, rightChoices = pfr, size = 8, multiple = F))
    })
    peak_prev_file("")
    peak_prev_name("")
  })
  
  observeEvent(peak_files, {
    if(length(peak_files()) == 0){
      output$SetChooser = renderUI({
        pfl = character()
        #add to active
        pfr = character()
        
        # print(pf)
        return(chooserInput(inputId = "ChooseIntervalSets", 
                            leftLabel = "Ready", rightLabel = "For Analysis", 
                            leftChoices = pfl, rightChoices = pfr, size = 8, multiple = F))
      })
    }
  })
  
  observeEvent(input$BtnCancelFile, {
    peak_prev_file("")
    peak_prev_name("")
  })
  
  observeEvent(input$FilesLoadData, {
    file_path = shinyFiles2load(input$FilesLoadData, roots_load)
    peak_prev_file(file_path)
    peak_prev_name(basename(file_path))
  })
  
  observeEvent(input$BtnUploadPeakfile, {
    peak_prev_file(input$BtnUploadPeakfile$datapath)
    peak_prev_name(input$BtnUploadPeakfile$name)
  })
  
 output$BtnDownloadResults = downloadHandler(
   filename = "intersectR.bed",
   content = function(con){
     write.table(gr_to_plot(), file = con, sep = "\t", col.names = F, row.names = F, quote = F)
   }
 )
  
  #whenever peak_prev_name changes (new file selected or file added/cancelled) update text box.
  observeEvent(peak_prev_name(), {
    new_val =  peak_prev_name()
    if(is.null(new_val)) new_val = ""
    # print(paste('update:', new_val))
    updateTextInput(session, "TxtFileName", value = new_val)
  })
  
  observeEvent(input$BtnDeleteSet, {
    sets = input$ChooseIntervalSets
    to_del = sets$selected
    
    sets$left = setdiff(sets$left, to_del)
    sets$right = setdiff(sets$right, to_del)
    output$SetChooser = renderUI({
      pfl = sets$left
      pfr = sets$right
      return(chooserInput(inputId = "ChooseIntervalSets", 
                          leftLabel = "Ready", rightLabel = "For Analysis", 
                          leftChoices = pfl, rightChoices = pfr, size = 8, multiple = F))
    })
    
  })
  observeEvent(input$BtnRenameSet, {
    sets = input$ChooseIntervalSets
    if(is.null(sets)) return(NULL)
    if(length(sets$selected) == 0) return(NULL)
    # print(sets)
    showModal(dataModal(sets = sets))
  })
  
  # When OK button is pressed, attempt to load the data set. If successful,
  # remove the modal. If not show another modal, but this time with a failure
  # message.
  observeEvent(input$BtnConfirmRename, {
    sets = input$ChooseIntervalSets
    old_name = sets$selected[1]
    new_name = input$TxtRename
    print(paste("rename", old_name, "-->", new_name))
    if(old_name == new_name){
      removeModal()
    }else if(any(c(unique(c(sets$left, sets$right)) == new_name))){
      showModal(dataModal(failed = TRUE))
    }else{
      tmp = peak_files()
      names(tmp)[names(tmp) == old_name] = new_name
      peak_files(tmp)
      
      
      sets$left[sets$left == old_name] = new_name
      sets$right[sets$right == old_name] = new_name
      output$SetChooser = renderUI({
        pfl = sets$left
        pfr = sets$right
        return(chooserInput(inputId = "ChooseIntervalSets", 
                            leftLabel = "Ready", rightLabel = "For Analysis", 
                            leftChoices = pfl, rightChoices = pfr, size = 8, multiple = F))
      })
      removeModal()
    }
  })
  
  output$NumericMergeExtensionOut = renderUI({
    if(!is.null(input$SliderMergeExtension)){
      num = input$SliderMergeExtension
    }else{
      num = 0
    }
    numericInput(inputId = "NumericMergeExtension", label = "", min = 0, max = Inf, value = num)
  })
  
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
  
  peak_dfr = reactive({
    if(is.null(peak_prev_file())) return(NULL)
    if(peak_prev_file() == "") return(NULL)
    load_peak_wValidation(peak_prev_file(), with_notes = T)
  })
  
  output$peaksHeader = DT::renderDataTable({
    if(is.null(peak_dfr())){
      m = matrix(0, ncol = length(peak_cn), nrow = 0)
      colnames(m) = peak_cn
      return(as.data.frame(m))
    }
    df = peak_dfr()
    # sdf = rbind(head(df), rep(".", ncol(df)), tail(df))
    DT::datatable(df, 
                  filter = list(position = "top", clear = TRUE, plain = F),
                  options = list(
                    pageLength = 5), rownames = F)
  })
  
  
  output$distPlot <- renderPlot({
    
    # generate bins based on input$bins from ui.R
    x    <- faithful[, 2]
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white')
    
  })
  
  observeEvent(input$BtnQuickFlat, handlerExpr = {
    peak_dirs = dir("/slipstream/galaxy/uploads/working/qc_framework/output_drugs_with_merged_inputs/", pattern = "MCF7_bza.+pooled", full.names = T)
    peak_files = sapply(peak_dirs[!grepl("_input_", peak_dirs)], function(d){dir(d, pattern = "peaks_passIDR", full.names = T)})
    names(peak_files) = sapply(strsplit(basename(peak_files), "_"), function(x)paste(x[1:3], collapse = "_"))
    peak_df = lapply(peak_files, load_peak_wValidation)
    for(i in 1:length(peak_df)){
      tmp = peak_files()
      tmp[[names(peak_df)[i]]] = peak_df[[i]]
      peak_files(tmp)
    }
    
    sets = input$ChooseIntervalSets
    output$SetChooser = renderUI({
      pfl = sets$left
      #add to active
      pfr = unique(c(sets$right, names(peak_df)))
      
      # print(pf)
      return(chooserInput(inputId = "ChooseIntervalSets", 
                          leftLabel = "Ready", rightLabel = "For Analysis", 
                          leftChoices = pfl, rightChoices = pfr, size = 8, multiple = F))
    })
    # peak_prev_file(NULL)
    # peak_prev_name(NULL)
  })
  
  gr_to_plot = reactiveVal(label = "gr_to_plot")
  
  #insulate other elements from needless updates when ChooseIntervalSets not really updated
  observeEvent(input$ChooseIntervalSets, {
    new_anlayze = input$ChooseIntervalSets$right
    names(new_anlayze) = new_anlayze
    was_analyze = to_analyze()
    if(length(was_analyze) != length(new_anlayze)){
      print("update anlaysis")
      to_analyze(new_anlayze)
    }else if(!all(was_analyze == new_anlayze)){
      print("update anlaysis")
      to_analyze(new_anlayze)
    }
    
  })
  
  to_analyze = reactiveVal(character())
  
  # observe(x = {
  #   to_analyze()
  #   input$StrategyRadio
  #   input$NumericMergeExtension
  #   gr_to_plot()
  #   peak_files()
  #   input$RadioPlotType
  #   showNotification("bwabwa")
  #   
  # })
  
  observeEvent(
    {
      #reactive dependencies
      to_analyze()
      input$StrategyRadio
      input$NumericMergeExtension
    }, {
      #prereqs
      if(length(to_analyze()) < 1) return(NULL)
      if(is.null(input$StrategyRadio)) return(NULL)
      if(is.null(input$NumericMergeExtension)) return(NULL)
      #body
      to_analyze_ = to_analyze()
      names(to_analyze_) = to_analyze_
      peak_df = peak_files()
      peak_df = peak_df[to_analyze_]
      peak_gr = lapply(peak_df, GRanges)
      print(names(peak_gr))
      showNotification("update gr_to_plot")
      gr_to_plot(intersectR(peak_gr, use_first = input$StrategyRadio == "serial", ext = input$NumericMergeExtension))
    })
  
  observeEvent({
    gr_to_plot()
    input$RadioPlotType
  }, {
    output$AnalysisPlot = renderPlot({
      isolate({
        print("try analysis render")
        
        if(is.null(input$NumericMergeExtension) || 
           is.null(input$ChooseIntervalSets)) return(NULL)
        if(length(input$ChooseIntervalSets$right) == 0) return(NULL)
        if(is.null(gr_to_plot())) return(NULL)
        
        gr = gr_to_plot()
        nam = names(peak_files())
        df = as.data.frame(elementMetadata(gr))
        tp = which(colnames(df) != "group")
        print("plot analysis render")
        showNotification("update plot", duration = .5)
        save(gr, df, tp, file = "last_plot.save")
        if(length(tp) < 2){
          plot(0:1, 0:1)
          text(.5, .5, "not enough sets for heatmap")
          return()
        }
       
        switch(input$RadioPlotType,
               "bars" = {
                 print("plot: bars")
                 hit_counts = colSums(df[,tp])
                 hit_counts_df = data.frame(count = hit_counts, 
                                            group = factor(names(hit_counts), levels = names(hit_counts)))
                 bp1<- ggplot(hit_counts_df, aes(x=group, y=count, fill=group))+
                   labs(x = "") +
                   geom_bar(width = 1, stat = "identity") +
                   guides(fill = "none") + 
                   theme_bw() +
                   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), 
                         panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
                 bp1 + annotate("text", y = hit_counts / 2, x = 1:length(hit_counts), label = hit_counts)
               },
               "pie" = {
                 print("plot: pie")
                 hit_counts = (colSums(df[,tp]))
                 hit_counts_df = data.frame(count = hit_counts, 
                                            group = factor(names(hit_counts), levels = names(hit_counts)))
                 
                 bp2<- ggplot(hit_counts_df, aes(x="", y=count, fill=(group)))+
                   labs(x = "") +
                   geom_bar(width = 1, stat = "identity") +
                   guides(x = "none") +
                   theme(axis.text.x = element_blank(), 
                         panel.background = element_blank(), 
                         axis.ticks = element_blank()) +
                   coord_polar("y", start = 0) + scale_y_reverse()
                 hc = rev(hit_counts) / sum(hit_counts)
                 hc = c(0, cumsum(hc)[-length(hc)]) + hc / 2
                 bp2 + annotate(geom = "text", y = hc * sum(hit_counts), x = 1.1, label = rev(hit_counts_df$count))
               },
               "venn" = {
                 print("plot: venn")
                 gg_color_hue <- function(n) {
                   hues = seq(15, 375, length = n + 1)
                   hcl(h = hues, l = 65, c = 100)[1:n]
                 }
                 source("jrb_vennDiagram.R")
                 jrb_venn(df[,tp], circle.col = gg_color_hue(length(tp)), cex = c(.8,1,1), bty = "n", names = "")
                 legend(x = "top", fill = gg_color_hue(length(tp)), legend = colnames(df)[tp], ncol = 2)
               },
               "euler" = {
                 print("plot: euler")
                 gg_venneuler(df[,tp])
               },
               "heatmap" = {  
                 print("plot: heatmap")
                 
                 # mat = mat + runif(length(mat), min = -.02, .02)
                 require(gplots)
                 if(length(tp) < 2){
                   plot(0:1, 0:1)
                   text(.5, .5, "not enough sets for heatmap")
                 }else{
                   # heatmap.2(mat[sample(nrow(mat), 2000),, drop = F]+.2, col = c("white", "black"), trace = "n")  
                   mat = ifelse(as.matrix(df[,tp]), 1, 0)
                   for(i in rev(1:ncol(mat))){
                     mat = mat[order(mat[,i]),]
                   }
                   hdt = as.data.table(cbind(mat, row = 1:nrow(mat)))
                   hdt = melt(hdt, id.vars = "row", variable.name = "groups")
                   hdt[, rnd := value > .5]
                   # png('tmp2.png', width = 1200, height = 1200)
                   require(png)
                   dmat = as.matrix(dcast(hdt, value.var = "value", formula = rev(row) ~ groups))[,1+1:length(tp)]
                   dmat = (dmat * -.95 + .95)
                   # amat = array(dmat, dim = c(nrow(dmat), ncol(dmat), 3) )
                   # colnames(dmat) = NULL
                   # rasterImage(dmat)
                   nr = 1000
                   nf = floor(nrow(dmat) / nr)
                   comp_mat = dmat[1:floor(nrow(dmat) / nf)*nf, rep(1:ncol(dmat), each = 20)]
                   writePNG(comp_mat, target = "tmp.png")
                   # writePNG(amat, target = "tmp.png")
                   require(grid)
                   png_grob = rasterGrob(readPNG("tmp.png"), width = unit(1, "npc"), height = unit(1, "npc"), interpolate = F)
                   # png_grob = rasterGrob(readPNG(system.file("img", "Rlogo.png", package="png")), width = unit(1, "npc"), height = unit(1, "npc"))
                   p = ggplot(hdt) + 
                     geom_tile(aes(x = groups, y = row, fill = NA, col = NULL)) +
                     scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "#EFEFEF")) +
                     scale_alpha(0) +
                     labs(fill = "binding", y = "", title = "Marks per region clustering") +
                     # theme_bw() +
                     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), 
                          # plot.margin = margin(r = .2, unit = "npc"),
                           panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
                           plot.background = element_blank(), panel.background = element_blank(), axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
                           axis.text = element_text(size = rel(1.5)), 
                           legend.key.size = unit(.12, units = "npc"),
                           # legend.key = element_rect(size = unit(.26, units = "npc")), 
                           legend.text = element_text(size = rel(2)),
                           legend.title = element_text(size = rel(3)),
                          
                           # legend.key.size = element_text(size = rel(2.5)),
                           axis.title = element_blank()) +
                     annotation_custom(png_grob, xmin = .5, xmax = ncol(dmat) + .5, ymin = .5, ymax = nrow(dmat)+.5)
                     # p + annotation_custom(rasterGrob(readPNG("tmp.png")))
                     # theme_set(base_size = 12)
                   print(p)
                   # print("done heatmap")
                   # dev.off()
                 }
                 
               })
      })
    })
  })
  
  
  observeEvent(input$FilesSaveResults, {
    if(is.null(input$FilesSaveResults)) return(NULL)
    print(input$FilesSaveResults)
    file_path = shinyFiles2save(input$FilesSaveResults, roots_output)
    print(file_path)
    if(!grepl(".bed$", file_path)){
      file_path = paste0(file_path, ".bed")
    }
    write.table(gr_to_plot(), file = file_path, row.names = F, col.names = T, quote = F)
  })
})
