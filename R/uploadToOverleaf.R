library(BiocWorkflowTools)
workflow_dir <- file.path(getwd(), 'workflow')
uploadToOverleaf(files = workflow_dir,
                 openInBrowser = TRUE,
                 forceNewProject = T)
