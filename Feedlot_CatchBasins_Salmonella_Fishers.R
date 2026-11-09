#Make 2x2 matrix with Salmonella culture status for CB and feces
salm_stat_mat <- matrix(c(3, 6, 2, 1), nrow = 2, byrow = TRUE,
              dimnames = list(c("Pen (Salmonella Positive)", "Pen (Salmonela Negative)"),
                              c("CB (Salmonella positive)", "CB (Salmonella Negative)")))

salm_stat_mat

#Run fisher's test
fisher.test(salm_stat_mat)
