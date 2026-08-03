box::use(
  shinyalert[
    shinyalert
  ]
)

#' @export
welcome_message <- function(app_name) {
  shinyalert(
    paste("Welcome to", app_name, " <br><small>
          KinGenie is being developed by the Sample Preparation and
          Characterisation (SPC)Facility at EMBL Hamburg.
          Please contact us if you have any questions,
          suggestions or feature requests.
          All comments are welcome.
          You can reach out us via spc@embl-hamburg.de or
          oburastero@embl-hamburg.de<br>"),
    imageWidth = 180, imageHeight = 180,
    closeOnClickOutside = FALSE, closeOnEsc = FALSE,
    confirmButtonText = "I accept", size = "m",
    showCancelButton = TRUE, cancelButtonText = "I decline", html = TRUE,
    confirmButtonCol = "#8bb8e8",
    callbackR = function(x) {
      if (!x) welcome_message(app_name)
    }
  )
}
#' @export
pop_up_warning <- function(string, size = "m") {
  shinyalert(
    text = string, type = "warning",
    closeOnEsc = TRUE,
    closeOnClickOutside = TRUE,
    html = TRUE, size = size
  )
}
#' @export
pop_up_info <- function(string, size = "m") {
  shinyalert(
    text = string, type = "info",
    closeOnEsc = TRUE,
    closeOnClickOutside = TRUE,
    html = TRUE, size = size
  )
}
#' @export
pop_up_success <- function(string, size = "m") {
  shinyalert(
    text = string, type = "success",
    closeOnEsc = TRUE, closeOnClickOutside = TRUE, html = TRUE, size = size
  )
}

# Print message to user regarding the simulation model
#' @export
print_model_message <- function(model_type_sim, model_selected_sim) {
  c1 <- model_type_sim == "surface"

  model_typ_text <- ifelse(c1, "a surface-based binding model
  (constant ligand concentration)",
    "an in solution-based binding model (ligand depletion)"
  )

  c2 <- model_selected_sim %in% c(
    "one_site_induced_fit",
    "one_site_conformational_selection"
  )

  extra_txt <- ifelse(!c1 && c2,
    "Relaxation rates are computed using Equations taken from Fabian Paul,
    Thomas R. Weikl, 2016.<br>", ""
  )

  extra_txt2 <- ifelse(c1 && model_selected_sim == "one_site_induced_fit",
    "In this model we assume that the induced and intermediate complex produce
    the same signal.<br>", ""
  )

  if (model_selected_sim == "one_site") {
    pop_up_info(
      paste0("The simulation will be performed using ", model_typ_text, "<br>
            The chemical equation is P + L ⇄ PL.
            Please wait until the simulation is completed."),
      size = "l"
    )
  }

  if (model_selected_sim == "heterogeneous_analyte") {
    pop_up_info(
      paste0("The simulation will be performed using ", model_typ_text, "<br>
            The chemical equations are
            L + P<sub>1</sub> ⇄ P<sub>1</sub>L
            (binding to the analyte P<sub>1</sub>) and
            L + P<sub>2</sub> ⇄ P<sub>2</sub>L
            (binding to the analyte P<sub>2</sub>). <br>
            The chemical reactions are: <br>
            <span style='display: inline-block; text-align: center;'>
            k<sub>on,1</sub><br>
            L + P<sub>1</sub> → P<sub>1</sub>L <br>
            k<sub>off,1</sub><br>
            P<sub>1</sub>L → L + P<sub>1</sub> <br>
            k<sub>on,2</sub><br>
            L + P<sub>2</sub> → P<sub>2</sub>L <br>
            k<sub>off,2</sub><br>
            P<sub>2</sub>L → L + P<sub>2</sub> <br>
            </span> <br>
            Please wait until the simulation is completed."),
      size = "l"
    )
  }


  if (model_selected_sim == "one_site_mtl") {
    pop_up_info(
      paste0('<p style="line-height: 2;">
            The simulation will be performed using ', model_typ_text, '<br>
            The chemical equations are
            L<sub>bulk</sub> ⇄ L<sub>surface</sub>
            (mass transport step) and
            L<sub>surface</sub> + P ⇄ PL (binding step). <br>
            The chemical reactions are: <br>
            <span style="display: inline-block; text-align: center;">
            k<sub>tr</sub><br>
            L<sub>bulk</sub> → L<sub>surface</sub> <br>
            k<sub>tr</sub><br>
            L<sub>surface</sub> → L<sub>bulk</sub> <br>
            k<sub>on</sub><br>
            L<sub>surface</sub> + P → PL <br>
            k<sub>off</sub><br>
            PL → L<sub>surface</sub> + P <br>
            </span> <br>
            ', extra_txt, "
            Please wait until the simulation is completed.
            </p>"),
      size = "l"
    )
  }

  if (model_selected_sim == "one_site_induced_fit") {
    pop_up_info(
      paste0('<p style="line-height: 2;">
            The simulation will be performed using ', model_typ_text, '<br>
            The chemical equations are
            P + L ⇄ PL<sub>int</sub> (binding step) and
            PL<sub>int</sub> ⇄ PL (induced fit step). <br>
            The chemical reactions are: <br>
            <span style="display: inline-block; text-align: center;">
            k<sub>on</sub><br>
            P + L → PL<sub>int</sub> <br>
            k<sub>off</sub><br>
            PL<sub>int</sub> → P + L <br>
            k<sub>c</sub><br>
            PL<sub>int</sub> → PL <br>
            k<sub>rev</sub><br>
            PL → PL<sub>int</sub> <br>
            </span> <br>
            ', extra_txt, ",
            ", extra_txt2, "
            Please wait until the simulation is completed.
            </p>"),
      size = "l"
    )
  }

  if (model_selected_sim == "one_site_conformational_selection") {
    pop_up_info(
      paste0('<p style="line-height: 2;">
            The simulation will be performed using ', model_typ_text, '<br>
            The chemical equations are
            P<sub>inactive</sub> ⇄ P<sub>active</sub>
            (conformational change step)
            and P<sub>active</sub> + L ⇄ PL (binding step). <br>
            The chemical reactions are: <br>
            <span style="display: inline-block; text-align: center;">
            k<sub>on</sub><br>
            P<sub>active</sub> + L → PL <br>
            k<sub>off</sub><br>
            PL → P<sub>active</sub> + L <br>
            k<sub>c</sub><br>
            P<sub>inactive</sub> → P<sub>active</sub> <br>
            k<sub>rev</sub><br>
            P<sub>active</sub> → P<sub>inactive</sub> <br>
            </span> <br>
            ', extra_txt, "
            Please wait until the simulation is completed.
            </p>"),
      size = "l"
    )
  }

  if (model_selected_sim == "ligand_has_two_sites") {
    pop_up_info(
      paste0('<p style="line-height: 2;">
            The simulation will be performed using ', model_typ_text, '<br>
            The chemical equations are
            P + L ⇄ PL (first binding step) and PL + L ⇄ PL<sub>2</sub>
            (second binding step). <br>
            The chemical reactions are: <br>
            <span style="display: inline-block; text-align: center;">
            2*k<sub>on</sub><br>
            P + L → PL (first binding step) <br>
            k<sub>off</sub><br>
            PL → P + L (first binding step) <br>
            sqrt(sigma) * k<sub>on</sub><br>
            PL → PL<sub>2</sub> (second binding step) <br>
            k<sub>off</sub> / sqrt(sigma)<br>
            PL<sub>2</sub> → PL (second binding step) <br>
            </span> <br>
            ', extra_txt, ",
            ", extra_txt2, "
            Please wait until the simulation is completed.
            </p>"),
      size = "l"
    )
  }

  NULL
}
