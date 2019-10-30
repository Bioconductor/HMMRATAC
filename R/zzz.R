.onLoad <- function(lib, pkg)
{
    rJava::.jpackage(pkg, jars = 'HMMRATAC.jar')
    rJava::.jaddClassPath(dir(file.path(getwd(), "inst/java"), full.names = TRUE))
    rJava::J("java.util.logging.LogManager")$getLogManager()$reset()
    invisible(NULL)
}
