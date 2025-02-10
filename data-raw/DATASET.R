
##Define example modules
example_fit<- fit_modules(example.voom)
example_wgna<- make_modules(example_fit, sft.value = 20)

example.mods<-
  example_wgna$mods%>%
  dplyr::arrange(module)

usethis::use_data(example.mods, overwrite = TRUE)
