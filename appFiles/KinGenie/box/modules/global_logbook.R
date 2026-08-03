box::use(
  shiny[
    moduleServer,
    reactive,
    reactiveVal
  ]
)

#' @export
logbookServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    logbook <- reactiveVal(character())

    append_record <- function(record,
                              include_time = FALSE,
                              add_empty_line = FALSE) {
      if (include_time) {
        record <- paste(
          format(Sys.time(), usetz = TRUE),
          record
        )
      }

      if (add_empty_line) {
        record <- c("", record)
      }

      logbook(c(logbook(), record))
    }

    clear <- function() {
      logbook(character())
    }

    list(
      logbook = reactive(logbook()),
      append = append_record,
      clear = clear
    )
  })
}
