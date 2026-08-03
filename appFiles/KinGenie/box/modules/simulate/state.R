box::use(
  shiny[
    moduleServer,
    reactiveValues
  ]
)

#' @export
simStateServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    reactiveValues(
      is_single_cycle_sim = FALSE,
      model_type_sim = "surface",
      model_selected_sim = "one_site",
      time_step = 0.5,
      init_lig_sim_surface = 0,
      init_lig_sim_solution = 0,
      protein_conc_sim = 0,
      protein_smax_sim = 0,
      numb_dil_sim = 0,
      lig_dil_factor_sim = 2,
      numb_dil_sim_prot = 0,
      prot_dil_factor_sim = 2,
      association_time = 300,
      dissociation_time = 600,
      total_time = 300,
      total_smax_sim = 5,
      kd_sim_1to1 = 0.5,
      koff_sim_1to1 = 0.01,
      ktr_sim_1to1 = 0.005,
      kon_sim_ligand_two_sites = 1e5,
      koff_sim_ligand_two_sites = 0.01,
      coop_sim_ligand_two_sites = 1,
      kd1_sim_hetero = 0.5,
      koff1_sim_hetero = 0.01,
      kd2_sim_hetero = 0.5,
      koff2_sim_hetero = 0.01,
      pop1_sim_hetero = 0.5,
      pop1_sim_smax = 0.5,
      kon_sim_1to1_adv = 1e5,
      koff_sim_1to1_adv = 0.01,
      kc_sim_1to1_adv = 0.005,
      krev_sim_1to1_adv = 0.005,
      signal_E = 1,
      signal_S = 1,
      signal_ES_simple = 1,
      signal_ES_int = 1,
      signal_ES = 1,
      pl_rmax_sim_two_sites = 1,
      pl2_rmax_sim_two_sites = 1
    )
  })
}
