calculate_evenness_metrics <- function(freqs) {
  # Number of alleles (S)
  S <- length(freqs)
  
  # Shannon-Wiener Index (H')
  H <- -sum(freqs * log(freqs))
  
  # Pielou's Evenness Index (J')
  J <- H / log(S)
  
  # Simpson's Diversity Index (D)
  D <- sum(freqs^2)
  
  # Simpson's Evenness Index (E1/D)
  E1_D <- (1 / D) / S
  
  # Shannon's Equitability Index (EH)
  EH <- H / log(S)
  
  # Berger-Parker Index
  p_max <- max(freqs)
  Berger_Parker <- 1 / p_max
  
  # Smith & Wilson’s Evenness Index (Evar)
  Evar <- 1 - (2 / pi) * atan(var(freqs) / mean(freqs))
  
  # Return a list of evenness metrics
  return(list(J = J, E1_D = E1_D, EH = EH, Berger_Parker = Berger_Parker, Evar = Evar))
}