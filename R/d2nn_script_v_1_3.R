d2nn_mono <- function (diretorio,verificado=F, processadores=T) {
  start_time_user <- Sys.time()
  
  suppressWarnings({
    #0.Lib----
    for (package in c('oro.dicom', 'oro.nifti','stringr','data.table','doParallel' ,'doSNOW','foreach')) {
      if (!require(package, character.only=T, quietly=T)) {
        install.packages(package)
        library(package, character.only=T)
      }
    }
    
    if(processadores==T){processadores<-detectCores()} else{processadores=processadores}
    
    #1.Definicao do diretorio padrao pelo usuário----
    #~/Projects/Neurocardio/GitHub/CONVERTER_DICOM_NIfTI/pasta_nao_compartilhada_somente_teste_BD/BB01
    cat("----------------------------------------------------------------------------------------------------")
    cat("\n\nBem vindo ao Conversor DICOM_2_NIfTI_Neurocardio_v1.3 [ D2NN_v1.3 ]\n\nAutor: Eric Wittlich  Data:14-07-2020  email:ericwittlich@gmail.com\n")
    dir_temp <-0
    if(verificado==F) {dir_ant = readline(cat("\n\nInstruções para uso do conversor:\n  1.Crie uma pasta num diretório conhecido, com o nome do DataSet (ex: C:/User/Fulano/Desktop/DataSet01). \n  2.Extraia todo o conteúdo do DataSet compactado para dentro da pasta criada no passo 1.\n  3.Escolha o método de entrada do endereço do DataSet a ser convertido:\n    > DIRETO: pressione ENTER e siga os comandos no 'Console'.\n    > ASSISTIDO: Utilize o painel 'File' (por default, o painel inferior direito do Rstudio) e siga os passos 'a' e 'b':\n      a) Localize o diretório através do ícone '...' na extrema direita do menu do painel, e abra-o.\n      b) Utilize o ícone da engrenagem ('More') no menu do painel e clique em 'Set as Working Directory'.\n  4.O processo será iniciado na sequência. Ao final será criada uma nova pasta, no mesmo local da pasta original, com a seguinte estrutura:\n      >Diretório original escolhido (ex: Desktop)\n         >Nome_do_DataSet_CLEANED(Pasta raiz dos novos arquivos)\n            >Nome_do_DataSet_NIfTI-files (Pasta onde estão os arquivos convertidos)\n            >AXIAL_T1/T2/FLAIR/ETC (Pastas com os arquivos DICOM originais de cada uma das sequências convertidas)\n\nATENÇÃO:\n > A conversão de cada slice leva, em média, 2 segundos. A montagem completa dos arquivos durará, em média, de 2 a 6 minutos - a depender da máquina e da resolução do arquivo.\n > O processo pode ser interrompido a qualquer instante utilizando o botão 'PARE' do 'Console' ou utilizando a tecla 'ESC'.\n > Caso ocorra algum erro durante o processo, ele será informado no 'Console'.\n > Caso deseje reiniciar o processo, por qualquer motivo, apague os arquivos gerados/duplicados e refaça as etapas desde o 'Passo 1'. \n\n----- Entre com o comando desejado do 'Passo 3' para continuar -----\n"))
    if(substr(dir_ant,0,5) == 'setwd'){dir_temp=(sub('.*"(.*)".*', "\\1", dir_ant))
    if(dir_temp!=0){diretorio=dir_temp}
    }
    }
    else{verificado==T}
    
    if(missing(diretorio)) {diretorio = readline(cat("\n>>> Defina a pasta do DataSet no 'Console'e pressione ENTER:"))} else {diretorio=diretorio}
    
    start_time <- Sys.time()
    cat("\nProcesso iniciado... Aguarde\n\n")
    
    #1.5Definicao do WD----
    
    
    setwd(diretorio)
    
    #2.Criacao do diretorio de trabalho, somente com arquivos DICOM----
    folder_name=gsub("^.*/","",diretorio)
    up_diretorio=str_remove(diretorio, folder_name)
    new_folder_name=paste(folder_name,"CLEANED",sep="_")
    new_diretorio=paste(up_diretorio,new_folder_name,sep="")
    to<-new_diretorio;from<-diretorio
    # Procura pelos arquivos DICOM no BBX
    files <- list.files(from, ".dcm" ,r = T)
    # criacao do BBX_CLEANED
    dir.create(to)
    # Cópia dos arquivos selecionados do BBX para o BBX_CLEANED
    file.copy(paste(from, files, sep = '/'), paste(to, files, sep = '/'))
    
    cat("+ OK - Criação do diretório de trabalho \n")
    
    #3.Criação da Matriz e separacao dos arquivos "axiais" [nome do arquivo,tipo de sequencia]----
    cat("--- Verificando consistência dos arquivos - Aguarde\n")
    
    Lista_Lv0<-data.frame(list.files(new_diretorio))
    Lista_Lv1 = NULL;fname=NULL;fhdr=NULL
    
    cl <- makeCluster(processadores)
    registerDoSNOW(cl)
    #progressbar
    
    pb <- txtProgressBar(max = (nrow(Lista_Lv0)), style = 3)
    progress <- function(n) setTxtProgressBar(pb,n)
    opts <-list(progress=progress)
    
    Lista_Lv1 <- foreach(i=1:(nrow(Lista_Lv0)), .combine = rbind, .options.snow = opts) %dopar% {
      fname[[i]] <- oro.dicom::readDICOM(Lista_Lv0[i,])
      fhdr[[i]] <- oro.dicom::extractHeader(fname[[i]]$hdr, string = "SeriesDescription", numeric = F, names = T, inSequence = T)
      axteste<-ifelse(grepl("axial",fhdr[[i]],ignore.case = T,perl = T,fixed=F),1,0)
      if (axteste==1) {
        Lista_Lv1[[i]] <- data.frame(Lista_Lv0[i,],fhdr[[i]])
      }
    }
    close(pb)
    stopCluster(cl)
    
    cat("\n+ OK - Verificação de consistência dos arquivos\n")
    
    #4.criacao das pastas AXIAL_T1, AXIAL_T2, AXIAL_FLAIR----
    T1dir<-paste(new_diretorio,"AXIAL_T1",sep="/")
    dir.create(T1dir)
    T2dir<-paste(new_diretorio,"AXIAL_T2",sep="/")
    dir.create(T2dir)
    FLAIRdir<-paste(new_diretorio,"AXIAL_FLAIR",sep="/")
    dir.create(FLAIRdir)
    
    #5.transfere os arquivos para as pastas AXIAL_T1, AXIAL_T2, AXIAL_FLAIR e limpa o diretorio----
    for (j in 1:(nrow(Lista_Lv1))) {
      T1_Lv1<-ifelse(grepl("T1",Lista_Lv1[[j,2]], ignore.case = T,perl = T, fixed = F),1,0)
      if (T1_Lv1==1) {
        file.copy(paste(new_diretorio,Lista_Lv1[[j,1]],sep = "/" ), paste(T1dir, Lista_Lv1[[j,1]], sep = "/"))
        file.remove(paste(new_diretorio,Lista_Lv1[[j,1]],sep = "/"))
      }
      else {
        next
      }
    }
    cat("+ OK - Diretório 'AXIAL_T1'\n")
    for (j in 1:(nrow(Lista_Lv1))) {
      T2_Lv1<-ifelse(grepl("T2",Lista_Lv1[[j,2]], ignore.case = T,perl = T, fixed = F),1,0)
      if (T2_Lv1==1) {
        file.copy(paste(new_diretorio,Lista_Lv1[[j,1]],sep = "/" ), paste(T2dir, Lista_Lv1[[j,1]], sep = "/"))
        file.remove(paste(new_diretorio,Lista_Lv1[[j,1]],sep = "/"))
      }
      else {
        next
      }
    }
    cat("+ OK - Diretório 'AXIAL_T2'\n")
    for (j in 1:(nrow(Lista_Lv1))) {
      FLAIR_Lv1<-ifelse(grepl("FLAIR",Lista_Lv1[[j,2]], ignore.case = T,perl = T, fixed = F),1,0)
      if (FLAIR_Lv1==1) {
        file.copy(paste(new_diretorio,Lista_Lv1[[j,1]],sep = "/" ), paste(FLAIRdir, Lista_Lv1[[j,1]], sep = "/"))
        file.remove(paste(new_diretorio,Lista_Lv1[[j,1]],sep = "/"))
      }
      else {
        next
      }
    }
    cat("+ OK - Diretório 'AXIAL_FLAIR'\n")
    lista_para_remover<-list.files(new_diretorio, ".dcm" ,r = F, include.dirs = F)
    file.remove(paste(new_diretorio,list.files(new_diretorio[4:length(list.files(lista_para_remover))]),sep = "/"))
    file.remove(diretorio)
    cat("+ OK - Limpeza de obsoletos \n")
    
    #6.Criacao dos arquivos NIfTI----
    
    cat("--- Início da conversão... Aguarde\n")
    
    niftiT1<-dicom2nifti(readDICOM(T1dir))
    cat("+ OK - T1_NIfTI\n")
    orthographic(niftiT1)
    cat("+ OK - T1_ortograf\n")
    
    niftiT2<-dicom2nifti(readDICOM(T2dir))
    cat("+ OK - T2_NIfTI\n")
    orthographic(niftiT2)
    cat("+ OK - T2_ortograf\n")
    
    niftiFLAIR<-dicom2nifti(readDICOM(FLAIRdir))
    cat("+ OK - FLAIR_NIfTI\n")
    orthographic(niftiFLAIR)
    cat("+ OK - FLAIR_ortograf\n")
    
    #7.Criacao da pasta nome_do_dataset_NIfTI_files e salvar arquivos NIfTI nesta pasta----
    niftidir<-paste(new_diretorio,paste(folder_name,"NIfTI_files",sep="_"),sep="/")
    dir.create(niftidir)
    writeNIfTI(niftiT1, paste(niftidir,paste(folder_name,"NIfTI_T1",sep="_"),sep="/"), onefile = TRUE, gzipped = TRUE, verbose = FALSE, warn = -1)
    cat("+ OK - T1_NIfTI armazenado\n")
    writeNIfTI(niftiT2, paste(niftidir,paste(folder_name,"NIfTI_T2",sep="_"),sep="/"), onefile = TRUE, gzipped = TRUE, verbose = FALSE, warn = -1)
    cat("+ OK - T2_NIfTI armazenado\n")
    writeNIfTI(niftiFLAIR, paste(niftidir,paste(folder_name,"NIfTI_FLAIR",sep="_"),sep="/"), onefile = TRUE, gzipped = TRUE, verbose = FALSE, warn = -1)
    cat("+ OK - FLAIR_NIfTI armazenado\n")
    
    #8.Apaga a pasta original dos arquivos descompactados e os arquivos----
    remfil <- list.files(diretorio,r = T, include.dirs = F, all.files = F)
    file.remove(remfil)
    file.remove(diretorio)
    
    
    #9.Messagem final----
    
    cat("\n+++++ Processo concluído +++++ \n\n")
    cat("Verifique os arquivos convertidos, em: ",new_diretorio,"\n")
    cat("\n")
    
  })
  cat("----------------------------------------------------------------------------------------------------\n")
  #cat("Tempos do conversor:\n")
  end_time <- Sys.time()
  
  cat("Number of Cores under load:",processadores,"\n")
  #useret = round(end_time - start_time_user, 2)
  useret = round(difftime(end_time,start_time_user,units = "secs"),0)
  #processet = round(end_time - start_time,2)
  processet = round(difftime(end_time,start_time,units = "secs"),0)
  cat("PROCESS elapsed time (s):",processet,"\n")
  cat("TOTAL elapsed time [user + process] (s):",useret,"\n")
  
}
