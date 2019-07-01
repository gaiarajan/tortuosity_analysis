# # Initialize empty vectors for vessel_name, resolution, tolerances, error_tot,
# error_ave, error_max, soam, dm, icm, tortuosity. For loop: from 0 to max
# number of vessels to be looked @ (5 for now)—> increment by 1 {if statement:
# !data points >20{continue} start=vessel_smooth.length*0.01; end =nrow things #
# save curvature vs torsion graph # run ivp, blah blah # go straight to 0.0001,
# save snapshot # add row # } # save csv
setwd("/Users/gayatrirajan/Desktop/externship/")
source("dat_file_reader.R")
source("single_vessel_plotter.R")
source("vessel_poly_fit.R")
source("vessel_spline_fit.R")
source("frenet_vectors.R")
source("tortuosity_metrics.R")
source("curve_torse_check_plotter.R")
source("subsample.R")
source("ivp_method_03.R")
library("rgl")
library("mgcv")
library("nat")
filename <-
  "MCAO day 7 microspheres-pecam NIH 1_vd_410x410x1000_561_sm_th_0.32.dat"
#the next bit defines zoom/viewing window params, which sets the viewing window. Just trying to minimize manual work.
vessels_slice <- dat_file_reader(dat_filename = filename)
vector1 = c(0.83157337, 0.09092994,-0.54792106, 0.00000000)
vector2 = c(-0.0197598, 0.9907266, 0.1344263, 0.0000000)
vector3 = c(0.5550634,-0.1009585, 0.8256586, 0.0000000)
vector4 = c(0, 0, 0, 1)
windowRect = c(0, 44, 512, 556)
sd = 100
zoom = 1
done = FALSE
userMatrix = array(c(vector1, vector2, vector3, vector4), dim = c(4, 4))
ave=c(1049.27,2048.53,1320.87,2742.05,1555.57,2272.35,2605.34,636.717,2425.41,2250.34,1276.33,1631.38,1326.61,1053.23,1703.47,718.066,1405.03,961.214,1042.61,1007.57,1389.43,1355.91,2261.99,2276.43,968.068,1019.94,1897.04,1681.59,987.617,1193.85,2008.73,2873.08,1873.7,2202.5,1841.31,1680.06,1672.96,1211.33,1611.61,905.044,2307.4,1718.65,1646.86,3069.94,1327.32,2530.55,867.395,880.122,1106.91,1483.59,1344.66,892.33,1258.18,1381.59,1916.61,1440.44,1297.66,1621.89,1551.82,1309.47,1388.35,656.858,1126.66,1502.41,1602.84,2144.38,878.322,1281.03,1042.18,1239.49,1364.68,2083.41,768.108,1397.85,1591.14,2093.73)
for (f in 1:76) {
  vessel <- vessels_slice[which(vessels_slice$ID == f), 1:3]
  if (length(vessel[[1]]) < 20) {
    next
  }
  vessel <- frenet_frame_calc(vessel_coords = vessel)
  error_bound = 50 * nrow(vessel[[1]])
  vessel_smooth <-
    vessel_spline_fit(
      vessel = vessel[[4]],
      number_samples = 20000,
      spline = "pspline",
      m = c(4, 2),
      subsample_density = sd,
      aic_slope = 15,
      plot = FALSE
    )
  while (((100 * nrow(vessel[[1]])) - nrow(vessel_smooth[[4]])) > error_bound &
         sd < 20000) {
    sd = sd + 185
    vessel_smooth <-
      vessel_spline_fit(
        vessel = vessel[[4]],
        number_samples = 20000,
        spline = "pspline",
        m = c(4, 2),
        subsample_density = sd,
        aic_slope = 15,
        plot = FALSE
      )
  }
  start <- nrow(vessel_smooth[[4]]) * 0.05
  end <- nrow(vessel_smooth[[4]])
  curve_torse_check_plotter(
    vessel_coords = vessel_smooth[[4]][start:end,],
    save_plot = TRUE,
    filename = paste("vessel_id_0",f,"_100X_sampling"),
    main = paste(
      "Curvature and Torsion vs Normalized Arclength \n Vessel ID",
      (f)
    ),
    col = "black",
    pch = 20
  )
  vessel_cur_tor_met <-
    curvature_torsion_calculator(vessel_coords = vessel_smooth[[4]][start:end,])
  vessel_ivp <-
    ivp_method(
      smth_vessel = vessel_smooth[[4]][start:end, ],
      tolerance1 = 0.0001,
      mkplot = TRUE,
      projection = TRUE,
      zoom = zoom,
      userMatrix = userMatrix,
      windowRect = windowRect
    )
  rgl.snapshot(filename = paste("vessel_id_0",f,"_100X_sampling_tol_0.0001.png"), fmt = "png")
  if (!done)
  {
    vessel_name <- paste("vessel_id_0", (f))
    resolution <- c("100X")
    tolerances <- c(0.0001)
    error_tot <- c(vessel_ivp[[3]])
    error_ave <- c(mean(vessel_ivp[[2]]))
    error_max <- c(max(vessel_ivp[[2]]))
    soam <-
      c(sum_of_all_angles_metric(vessel_coords = vessel_smooth[[4]])[[1]][[1]])
    dm <-
      c(distance_metric(vessel_coords = vessel_smooth[[4]])[2][[1]])
    icm <-
      c(inflection_count_metric(vessel_coords = vessel_smooth[[4]]))
    ave_radius <- ave[[f]]/1000
    tortuosity <-
      data.frame(
        "total_torsion" = c(vessel_cur_tor_met$TT),
        "max_torsion" = c(vessel_cur_tor_met$MT),
        "total_curvature" = c(vessel_cur_tor_met$TC),
        "max_curvature" = c(vessel_cur_tor_met$MC),
        "average_curvature" = c(vessel_cur_tor_met$AC),
        "average_torsion" = c(vessel_cur_tor_met$AT),
        "total_combined_curvature" = c(vessel_cur_tor_met$TCC),
        "average_combined_curvature" = c(vessel_cur_tor_met$ACC),
        "total_squared_curvature" = c(vessel_cur_tor_met$TCsq),
        "total_squared_torsion" = c(vessel_cur_tor_met$TTsq)
      )
    done=TRUE
  }
  else
  {
    vessel_name <- c(vessel_name, paste("vessel_id_", (f)))
    resolution <- c(resolution, "100X")
    tolerances <- c(tolerances, 0.0001)
    error_tot <- c(error_tot, vessel_ivp[[3]])
    error_ave <- c(error_ave, mean(vessel_ivp[[2]]))
    ave_radius <- c(ave_radius,ave[[f]]/1000)
    error_max <- c(error_max, max(vessel_ivp[[2]]))
    soam <-
      c(soam,
        sum_of_all_angles_metric(vessel_coords = vessel_smooth[[4]][start:end, ])[[1]][[1]])
    dm <-
      c(dm, distance_metric(vessel_coords = vessel_smooth[[4]][start:end, ])[2][[1]])
    icm <-
      c(icm,
        inflection_count_metric(vessel_coords = vessel_smooth[[4]][start:end, ]))
    tortuosity <-
      rbind(
        tortuosity,
        data.frame(
          "total_torsion" = c(vessel_cur_tor_met$TT),
          "max_torsion" = c(vessel_cur_tor_met$MT),
          "total_curvature" = c(vessel_cur_tor_met$TC),
          "max_curvature" = c(vessel_cur_tor_met$MC),
          "average_curvature" = c(vessel_cur_tor_met$AC),
          "average_torsion" = c(vessel_cur_tor_met$AT),
          "total_combined_curvature" = c(vessel_cur_tor_met$TCC),
          "average_combined_curvature" = c(vessel_cur_tor_met$ACC),
          "total_squared_curvature" = c(vessel_cur_tor_met$TCsq),
          "total_squared_torsion" = c(vessel_cur_tor_met$TTsq)
        )
      )
  }
}
summary_stats <- data.frame("vessel" = vessel_name, "resolution" = resolution, "tolerances" = tolerances, "error_tot" = error_tot, "error_ave" = error_ave, "error_max" = error_max, "curvature_total" = tortuosity$total_curvature, "curvature_max" = tortuosity$max_curvature, "curvature_average" = tortuosity$average_curvature, "torsion_total" = tortuosity$total_torsion, "torsion_max" = tortuosity$max_torsion, "torsion_average" = tortuosity$average_torsion, "total_combined_curvature" = tortuosity$total_combined_curvature, "average_combined_curvature" = tortuosity$average_combined_curvature, "total_squaredcurvature" = tortuosity$total_squared_curvature, "total_squared_torsion" = tortuosity$total_squared_torsion, "SOAM" = soam, "DM" = dm, "ICM" = icm, "min_radius"=ave_radius)

write.csv(x = summary_stats, file = "torsion_versus_resolution_tolerance_summary_stats.csv")
