# Implementation of the One Sample Test
#
# This file is part of the tTestLFHK jamoovi package.
#
# Copyright 2022--2026 Gisbert W. Selke <gwselke@github.com>
#
# This implementation is NOT optimized.
# * It could make use of jamovi's clearWith feature (but it is a bit tricky to get this right for all cases).
# * Shapiro-Wilks test could be cached (per variable), because it is independent of all other settings.
# * t-test/Wilcoxon test results could be cached (per variable, confLevel, mean, and alternative).
# But all the tests are fast, so optimisation probably does not pay off.
#

ttestOSClass <- if ( requireNamespace('jmvcore', quietly=TRUE) ) R6::R6Class(
  "ttestOSClass",
  inherit = ttestOSBase,
  private = list(

    .init = function( ) {
      # Do some initialisations on the result tables that could not be done statically via ttestos.a.yaml

      isParametric <- self$options$testType == 'parametric'
      tableT       <- if (isParametric) self$results$ttestParOS else self$results$ttestNonparOS
      tableS       <- self$results$normalityOS

      if (isParametric) {
        confLevel      <- self$options$confLevel
        tableT$addColumn( name        = 'mean_mu0',
                          title       = 'Mean - μ<sub>0</sub>',
                          type        = 'number'
                        )
        superTitleText <- paste0( 100*confLevel, '% Confidence Interval' )                        
        tableT$addColumn( name        = 'cilower',
                          title       = 'Lower',
                          superTitle  = superTitleText,
                          type        = 'number'
                        )
        tableT$addColumn( name        = 'ciupper',
                          title       = 'Upper',
                          superTitle  = superTitleText,
                          type        = 'number'
                        )
        tableT$setNote( 'testvalue', 'Test value = hypothesized population mean (μ<sub>0</sub>)', init=TRUE )
        tableT$setNote( 'nullhyp', paste( 'H<sub>0</sub>: μ =', self$options$mean, '(i.e. the population mean μ is equal to the hypothesized value μ<sub>0</sub>)' ), init=TRUE )
        tableT$setNote( 'althyp', 
                        paste( 'H<sub>A</sub>:  μ',
                               if ( self$options$alternative=='two.sided' ) '≠' else if ( self$options$alternative=='less' ) '<' else '>',
                               self$options$mean
                             ), 
                        init=TRUE  
                      )
        tableT$setNote( 'significance', 'If the p-value is smaller than the significance level, we reject H<sub>0</sub>.', init=TRUE )
      } else {
        tableT$setNote( 'testvalue', 'Test value = hypothesized centre of symmetry of the distribution', init=TRUE )
        tableT$setNote( 'nullhyp', 
                        paste0( 'H<sub>0</sub>: The centre of symmetry of the distribution is equal to ', self$options$mean, '.' ), 
                        init=TRUE 
                      )
      }
      
      tableS$setNote( 'normality', 'H<sub>0</sub>: the data come from a normal (i.e. Gaussian) distribution.', init=TRUE )
      tableS$setNote( 'normality_sig', 'If the p-value is smaller than the significance level, we reject H<sub>0</sub>.', init=TRUE )

      self$results$ttestParOS$setVisible(isParametric)
      self$results$ttestNonparOS$setVisible(!isParametric)
      tableS$setVisible( self$options$normalityTest )
    }, ## end .init

    .run = function( ) {

      # `self$data`    contains the data
      # `self$options` contains the options
      # `self$results` contains the results object (to populate)

      isParametric <- self$options$testType == 'parametric'
      tableT       <- if (isParametric) self$results$ttestParOS else self$results$ttestNonparOS
      tableS       <- self$results$normalityOS

      warn_nn_ct   <- 0
      # note_sig_ct  <- 0

      for ( thisx in self$options$xs ) {
        # for all chosen variables, calculate the appropriate test, including Shapiro-Wilks, and add a line to the result table

        x <- jmvcore::toNumeric( self$data[[thisx]] )
        if ( is.factor(x) ) jmvcore::reject( paste( 'Cannot run test on grouping variable', thisx ) )

        if (isParametric) {
          # Parametric case: Student's t
          results <- t.test( x           = x,
                             alternative = self$options$alternative,
                             mu          = self$options$mean,
                             conf.level  = self$options$confLevel,
                             paired      = FALSE
                           )

          df        <- results$parameter
          sided     <- if ( results$alternative=='two.sided' ) 2 else 1
          alpha     <- 1 - attr( results$conf.int, 'conf.level')
          critval   <- qt( 1 - alpha/sided, df )
          if ( results$alternative=='less' ) critval = -critval

          # Add row to main result table
          tableT$setRow( rowKey = thisx,
                         values = list(
                                        t        = results$statistic,
                                        testName = "Student's t",
                                        crit     = critval,
                                        df       = df,
                                        p        = results$p.value,
                                        mean     = results$estimate,
                                        mean_mu0 = results$estimate - self$options$mean,
                                        cilower  = results$conf.int[1],
                                        ciupper  = results$conf.int[2]
                                      )
                       )
        } else {
          # Non-parametric case: Wilcoxon
          results <- wilcox.test(
                                  x           = x,
                                  alternative = self$options$alternative,
                                  mu          = self$options$mean,
                                  digits.rank = 7,
                                  conf.level  = self$options$confLevel,
                                  conf.int    = TRUE,
                                  paired      = FALSE
                                )

          # Add row to main result table
          tableT$setRow( rowKey = thisx,
                         values = list(
                                        t        = results$statistic,
                                        testName = "Wilcoxon W",
                                        p        = results$p.value
                                      )
                       )
        } ## if (isParametric)

        # Test for normality; maybe add a warning marker to this variable:
        resultsS <- shapiro.test(x)
        if (isParametric) {
          if ( resultsS$p.value < 1-self$options$confLevel ) {
            tableT$addSymbol( rowKey=thisx, col=1, symbol='‡' )
            warn_nn_ct <- warn_nn_ct + 1
          } else {
            # If normality is not violated, possibly add a marker for significant result:
            # if ( results$p.value < 1-self$options$confLevel ) {
            #   tableT$addSymbol( rowKey=thisx, col=1, symbol='*' )
            #   note_sig_ct <- note_sig_ct + 1
            # }
          }
        }


        if ( self$options$normalityTest ) {
          # If user asked for display of normality test, add a line to normality test table:
          tableS$setRow( rowKey = thisx,
                         values = list(
                                        w = resultsS$statistic,
                                        p = resultsS$p.value
                                      )
                       )
        }

      }

      if ( warn_nn_ct > 0  ) tableT$setNote( 'normalWarning',
                                             paste( 'For variables marked ‡ the assumption of normality may be violated.',
                                                    'Recommendation: use a nonparametric test (One Sample Wilcoxon Test).'
                                                  ),
                                             init=FALSE
                                           )
      # if ( note_sig_ct > 0 ) tableT$setNote( 'significantNote1',
      #                                        paste( 'For variables marked * the calculated p-value is smaller than the chosen significance level.',
      #                                               'Therefore we reject the null hypothesis for each variable marked with *.'
      #                                             ),
      #                                        init=FALSE
      #                                      )

    } ## end .run

  )
)
