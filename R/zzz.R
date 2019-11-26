.onLoad <- function(lib, pkg)
{
    rJava::.jpackage(pkg, jars = 'HMMRATAC.jar')
    path <- list.files(system.file(package = pkg, "inst", "java"), full.names = TRUE)
    rJava::.jaddClassPath(path)
    rJava::J("java.util.logging.LogManager")$getLogManager()$reset()
    invisible(NULL)
}
