
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
source('server_setup.R')

shinyServer(function(input, output, session) {
  #ShinyFile initialization
  shinyFileChoose(input, 'FilesLoadData', roots= roots_load, filetypes=c("narrowPeak", "broadPeak", "bed", "txt"))
  shinyFileSave(input, 'FilesSaveResults', roots= roots_output, filetypes=c("bed"))
  
  #Set Preview reactives
  PreviewSet_Filepath = reactiveVal(value = "", label = "PreviewSet_Filepath")
  PreviewSet_Name = reactiveVal(value = "", label = "PreviewSet_Name")
  PreviewSet_DataFrame = reactive({
    if(is.null(PreviewSet_Filepath())) return(NULL)
    if(PreviewSet_Filepath() == "") return(NULL)
    load_peak_wValidation(PreviewSet_Filepath(), with_notes = T)
  })
  
  #Set Organization of Loaded reactives
  SetsLoaded_Selected = reactiveVal(value = character(), label = "SetsLoaded_Selected")
  SetsLoaded_DataFrames = reactiveVal(value = list(), label = "SetsLoaded_DataFrames")
  SetsLoaded_FilePaths = reactiveVal(value = list(), label = "SetsLoaded_FilePaths")
  
  
  observeEvent(input$BtnAddFile, {
    if(is.null(PreviewSet_DataFrame())) return(NULL)
    tmp = SetsLoaded_DataFrames()
    tmp2 = SetsLoaded_FilePaths()
    name_in_list = input$TxtFileName
    if(any(names(tmp) == name_in_list)){
      n = sum(grepl(name_in_list, names(tmp), fixed = T))
      name_in_list = paste0(name_in_list, "(", n, ")")
    }
    tmp[[name_in_list]] = PreviewSet_DataFrame()[input$peaksHeader_rows_all,]
    tmp2[[name_in_list]] = PreviewSet_Filepath()
    SetsLoaded_DataFrames(tmp)
    SetsLoaded_FilePaths(tmp2)
    sets = input$ChooseIntervalSets
    output$SetChooser = renderUI({
      pfl = sets$left
      #add to active
      pfr = c(sets$right, name_in_list)
      
      # print(pf)
      return(chooserInput(inputId = "ChooseIntervalSets", 
                          leftLabel = "Ready", rightLabel = "For Analysis", 
                          leftChoices = pfl, rightChoices = pfr, size = 8, multiple = F))
    })
    PreviewSet_Filepath("")
    PreviewSet_Name("")
  })
  
  observeEvent(SetsLoaded_DataFrames(), {
    if(length(SetsLoaded_DataFrames()) == 0){
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
    PreviewSet_Filepath("")
    PreviewSet_Name("")
  })
  
  observeEvent(input$FilesLoadData, {
    file_path = shinyFiles2load(input$FilesLoadData, roots_load)
    PreviewSet_Filepath(file_path)
    PreviewSet_Name(basename(file_path))
  })
  
  observeEvent(input$BtnUploadPeakfile, {
    PreviewSet_Filepath(input$BtnUploadPeakfile$datapath)
    PreviewSet_Name(input$BtnUploadPeakfile$name)
  })
  
  output$DownloadResults = downloadHandler(
    filename = "intersectR.bed",
    content = function(con){
      write.table(gr_to_plot(), file = con, sep = "\t", col.names = F, row.names = F, quote = F)
    }
  )
  output$DownloadPlot = downloadHandler(
    filename = "intersectR.pdf",
    content = function(con){
      
      pdf(con)
      if(is.null(last_plot())){
        plot(0:1, 0:1)
        text(.5, .5, "not enough data to plot")
      }else{
        print(last_plot())
      }
      dev.off()
    }
  )
  
  ##can't disable DL buttons for some reason
  # observe({
  #   if(is.null(last_plot())){
  #     showNotification("disable plot DL")
  #     disable(input$DownloadPlot)
  #   }
  #   if(!is.null(last_plot())){
  #     showNotification("enable plot DL")
  #     enable(input$DownloadPlot)
  #   }
  # })
  
  #whenever PreviewSet_Name changes (new file selected or file added/cancelled) update text box.
  observeEvent(PreviewSet_Name(), {
    new_val =  PreviewSet_Name()
    if(is.null(new_val)) new_val = ""
    # print(paste('update:', new_val))
    updateTextInput(session, "TxtFileName", value = new_val)
  })
  
  #keep setsloaded and such updatd
  observeEvent(
    eventExpr = input$ChooseIntervalSets, 
    handlerExpr = {
      sets = input$ChooseIntervalSets
      active = sets$right
      selected = paste("selected:", sets$selected)
      file_paths =  SetsLoaded_FilePaths()
      data_frames = SetsLoaded_DataFrames()
      if(length(active) > 0){
        active_list = sapply(1:length(active), function(i){
          paste0("\t", i, ") ", active[i], ": ", nrow(data_frames[[active[i]]]))
        })
      }else{
        active_list = character()
      }
      str = paste(sep = "\n",
                  "---Selection Debug Info---",
                  selected, 
                  "active sets:",
                  paste(active_list, collapse = "\n"))
      # showNotification(
      output$DebugTxt = renderText({
        return(str)
      })
      # )
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
  observeEvent(input$BtnFilterSet, {
    if(length(input$ChooseIntervalSets$selected) != 1){
      showNotification("No data selected.", type = "error")
      return()
    }
    showModal(filterModal())
  })
  observeEvent(input$BtnConfirmFilter, {
    showNotification("Confirm filter.", type = "message")
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
      #if no change, that's fine
      removeModal()
    }else if(any(c(unique(c(sets$left, sets$right)) == new_name))){
      #if duplicated name, go to failure
      showModal(dataModal(sets, failed = TRUE))
    }else{
      #do the renaming and update stuff
      tmp = SetsLoaded_DataFrames()
      names(tmp)[names(tmp) == old_name] = new_name
      SetsLoaded_DataFrames(tmp)
      
      
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
  
  
  

  
  output$DTPeaksHeader = DT::renderDataTable({
    if(is.null(PreviewSet_DataFrame())){
      m = matrix(0, ncol = length(peak_cn), nrow = 0)
      colnames(m) = peak_cn
      return(as.data.frame(m))
    }
    df = PreviewSet_DataFrame()
    # sdf = rbind(head(df), rep(".", ncol(df)), tail(df))
    DT::datatable(df, 
                  filter = list(position = "top", clear = TRUE, plain = F),
                  options = list(
                    pageLength = 5), rownames = F)
  })
  
  output$DTPeaksFilter = DT::renderDataTable({
    
    df = SetsLoaded_DataFrames()[[input$ChooseIntervalSets$selected]]
    # sdf = rbind(head(df), rep(".", ncol(df)), tail(df))
    DT::datatable(df, 
                  filter = list(position = "top", clear = TRUE, plain = F),
                  options = list(
                    scrollX = T,
                    pageLength = 10), rownames = F)
  })
  
  observeEvent(input$BtnQuickFlat, handlerExpr = {
    peak_dirs = dir("/slipstream/galaxy/uploads/working/qc_framework/output_drugs_with_merged_inputs/", pattern = "MCF7_bza.+pooled", full.names = T)
    peak_files = sapply(peak_dirs[!grepl("_input_", peak_dirs)], function(d){dir(d, pattern = "peaks_passIDR", full.names = T)})
    names(peak_files) = sapply(strsplit(basename(peak_files), "_"), function(x)paste(x[1:3], collapse = "_"))
    peak_df = lapply(peak_files, load_peak_wValidation)
    for(i in 1:length(peak_df)){
      tmp = SetsLoaded_DataFrames()
      tmp[[names(peak_df)[i]]] = peak_df[[i]]
      SetsLoaded_DataFrames(tmp)
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
  #   SetsLoaded_DataFrames()
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
      peak_df = SetsLoaded_DataFrames()
      peak_df = peak_df[to_analyze_]
      peak_gr = lapply(peak_df, GRanges)
      print(names(peak_gr))
      # showNotification("update gr_to_plot")
      gr_to_plot(intersectR(peak_gr, use_first = input$StrategyRadio == "serial", ext = input$NumericMergeExtension))
    })
  
  last_plot = reactiveVal(NULL)
  
  output$AnalysisPlot = renderPlot({
    p = last_plot()
    if(is.null(p)){
      plot(0:1, 0:1)
      text(.5, .5, "waiting on data to plot")
    }else{
      p
    }
  })
  
  # observe(
  #   {
  #     last_plot()
  #     showNotification("last plot changed")
  #   }
  # )
  
  observeEvent({
    gr_to_plot()
    input$RadioPlotType
  }, {
    # showNotification("try analysis plot render")
    if(is.null(input$NumericMergeExtension) || 
       is.null(input$ChooseIntervalSets)) return(NULL)
    if(length(input$ChooseIntervalSets$right) == 0) return(NULL)
    if(is.null(gr_to_plot())) return(NULL)
    
    gr = gr_to_plot()
    nam = names(SetsLoaded_DataFrames())
    df = as.data.frame(elementMetadata(gr))
    tp = which(colnames(df) != "group")
    print("plot analysis render")
    # showNotification("update plot", duration = .5)
    save(gr, df, tp, file = "last_plot.save")
    if(length(tp) < 2){
      plot(0:1, 0:1)
      text(.5, .5, "not enough sets for heatmap")
      p = NULL
    }
    
    switch(input$RadioPlotType,
           "bars" = {
             print("plot: bars")
             hit_counts = colSums(df[,tp, drop = F])
             hit_counts_df = data.frame(count = hit_counts, 
                                        group = factor(names(hit_counts), levels = names(hit_counts)))
             bp1<- ggplot(hit_counts_df, aes(x=group, y=count, fill=group))+
               labs(x = "") +
               geom_bar(width = 1, stat = "identity") +
               guides(fill = "none") + 
               theme_bw() +
               theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), 
                     panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
             bp1 = bp1 + annotate("text", y = hit_counts / 2, x = 1:length(hit_counts), label = hit_counts)
             p = bp1
           },
           "pie" = {
             print("plot: pie")
             hit_counts = (colSums(df[,tp, drop = F]))
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
             bp2 = bp2 + annotate(geom = "text", y = hc * sum(hit_counts), x = 1.1, label = rev(hit_counts_df$count))
             p = bp2
           },
           "venn" = {
             print("plot: venn")
             gg_cols = gg_color_hue(length(tp))
             if(length(tp) > 3){
               showNotification("venn currently limited to 3 sets!", type = "warning")
               tp = tp[1:3]
             }
             source("jrb_vennDiagram.R")
             p = gg_venn(df[,tp], labels_size = 8, counts_size = 8, circle.col = gg_cols[1:length(tp)])
             # v = jrb_venn(df[,tp], circle.col = gg_color_hue(length(tp)), cex = c(.8,1,1), bty = "n", names = "")
             # legend(x = "top", fill = gg_color_hue(length(tp)), legend = colnames(df)[tp], ncol = 2)
           },
           "euler" = {
             print("plot: euler")
             if(length(tp) > 1){
               p = gg_venneuler(memb = df[,tp, drop = F])
             }else{
               showNotification("euler plot for 1 group doesn't make sense.  here's a venn.", type = "warning")
               p = gg_venn(df[,tp], labels_size = 8, counts_size = 8, circle.col = gg_color_hue(length(tp)))
             }
           },
           "heatmap" = {  
             print("plot: heatmap")
             
             # mat = mat + runif(length(mat), min = -.02, .02)
             require(gplots)
             if(length(tp) < 1){
               plot(0:1, 0:1)
               text(.5, .5, "not enough sets for heatmap")
               p = NULL
             }else{
               mat = ifelse(as.matrix(df[,tp, drop = F]), 1, 0)
               for(i in rev(1:ncol(mat))){
                 mat = mat[order(mat[,i]), , drop = F]
               }
               hdt = as.data.table(cbind(mat, row = 1:nrow(mat)))
               hdt = melt(hdt, id.vars = "row", variable.name = "groups")
               hdt[, rnd := value > .5]
               require(png)
               dmat = as.matrix(dcast(hdt, value.var = "value", formula = rev(row) ~ groups))[,1+1:length(tp), drop = F]
               dmat = (dmat * -.95 + .95)
               nr = 1000
               nf = floor(nrow(dmat) / nr)
               comp_mat = dmat[1:floor(nrow(dmat) / nf)*nf, rep(1:ncol(dmat), each = 20)]
               writePNG(comp_mat, target = "tmp.png")
               require(grid)
               png_grob = rasterGrob(readPNG("tmp.png"), width = unit(1, "npc"), height = unit(1, "npc"), interpolate = F)
               p = ggplot(hdt) + 
                 geom_tile(aes(x = groups, y = row, fill = NA, col = NULL)) +
                 scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "#EFEFEF")) +
                 scale_alpha(0) +
                 labs(fill = "binding", y = "", title = "Marks per region clustering") +
                 theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), 
                       panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
                       plot.background = element_blank(), panel.background = element_blank(), axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
                       axis.text = element_text(size = rel(1.5)), 
                       legend.key.size = unit(.12, units = "npc"),
                       legend.text = element_text(size = rel(2)),
                       legend.title = element_text(size = rel(3)),
                       axis.title = element_blank()) +
                 annotation_custom(png_grob, xmin = .5, xmax = ncol(dmat) + .5, ymin = .5, ymax = nrow(dmat)+.5)
             }
           })
    last_plot(p)
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
