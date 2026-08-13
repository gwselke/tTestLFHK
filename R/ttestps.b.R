# Implementation of the paired sample test
#
# This file is part of the zTestLFHK jamoovi package.
#
# Copyright 2022--2026 Gisbert W. Selke <gwselke@github.com>
#
# This implementation is NOT optimized.
# * It could make use of jamovi's clearWith feature (but it is a bit tricky to get this right for all cases).
# * Shapiro-Wilks test could be cached (per variable), because it is independent of all other settings.
# * t-test results could be cached (per variable, confLevel, mean, and alternative).
# But all the tests are fast, so optimisation probably does not pay off.
#

ttestPSClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
  "ttestPSClass",
  inherit = ttestPSBase,
  private = list(

    .init = function( ) {
      # Do some initialisations on the result tables that could not be done statically via ttestps.a.yaml

      isParametric <- self$options$testType == 'parametric'
      tableT       <- if (isParametric) self$results$ttestParPS else self$results$ttestNonparPS
      tableS       <- self$results$normalityPS

      if (isParametric) {
        confLevel      <- self$options$confLevel
        superTitleText <- paste0( 100*confLevel, '% Confidence Interval' )

        tableT$addColumn( name       = 'cilower',
                          title      = 'Lower',
                          superTitle = superTitleText,
                          type       = 'number'
                        )
        tableT$addColumn( name       = 'ciupper',
                          title      = 'Upper',
                          superTitle = superTitleText,
                          type       = 'number'
                        )
                        
        tableT$setNote( 'nullhyp', 'H<sub>0</sub>: μ<sub>diff</sub> = 0 (i.e., the population mean of the differences is 0', init=TRUE )
        descr <- paste( 'H<sub>A</sub>: μ<sub>diff</sub>',
                        if ( self$options$alternative=='two.sided' ) '≠' else if ( self$options$alternative=='less' ) '<' else '>',
                        '0'
                      )
        tableT$setNote( 'althyp', descr, init=TRUE )
        tableT$setNote( 'significance', 'If the p-value is smaller than the significance level, we reject H<sub>0</sub>.', init=TRUE )
      }

      tableS$setNote( 'normality', 'H<sub>0</sub>: the distribution of the pairwise differences is normal (i.e. Gaussian).', init=TRUE )
      tableS$setNote( 'normality_sig', 'If the p-value is smaller than the significance level, we reject H<sub>0</sub>.', init=TRUE )

      self$results$ttestParPS$setVisible(isParametric)
      self$results$ttestNonparPS$setVisible(!isParametric)
      tableS$setVisible( self$options$normalityTest )

    }, ## end .init

    .run = function( ) {

      # `self$data`    contains the data
      # `self$options` contains the options
      # `self$results` contains the results object (to populate)

      # if ( is.null(self$options$x) || is.null(self$options$y) ) return( )

      isParametric <- self$options$testType == 'parametric'
      tableT       <- if (isParametric) self$results$ttestParPS else  self$results$ttestNonparPS
      tableS       <- self$results$normalityPS

      warn_nn_ct   <- 0
      # note_sig_ct  <- 0

      for ( key in tableT$rowKeys ) {

        xname <- key[[1]]
        yname <- key[[2]]
        if ( is.null(xname) || is.null(yname) ) next( )

        x <- jmvcore::toNumeric( self$data[[xname]] )
        y <- jmvcore::toNumeric( self$data[[yname]] )
        if ( is.factor(x) ) jmvcore::reject( paste( 'Cannot run test on grouping variable', xname ) )
        if ( is.factor(y) ) jmvcore::reject( paste( 'Cannot run test on grouping variable', yname ) )

        if (isParametric) {

          # Parametric case: Student's t test
          results <- t.test(
                             x           = x,
                             y           = y,
                             alternative = self$options$alternative,
                             mu          = 0,
                             conf.level  = self$options$confLevel,
                             paired      = TRUE
                           )
          df      <- results$parameter
          sided   <- if ( results$alternative=='two.sided' ) 2 else 1
          alpha   <- 1 - attr( results$conf.int, 'conf.level')
          critval <- qt( 1 - alpha/sided, df )
          if ( results$alternative=='less' ) critval = -critval

          # Add row to main result table
          tableT$setRow( rowKey = key,
                         values = list(
                                        x        = xname,
                                        y        = yname,
                                        testName = "Student's t",
                                        t        = results$statistic,
                                        crit     = critval,
                                        df       = df,
                                        p        = results$p.value,
                                        mean     = results$estimate,
                                        cilower  = results$conf.int[1],
                                        ciupper  = results$conf.int[2]
                                      )
                       )
        } else {

          # Nonparametric case: Wilcoxon W test
          results <- wilcox.test(
                                  x           = x,
                                  y           = y,
                                  alternative = self$options$alternative,
                                  mu          = 0,
                                  conf.level  = self$options$confLevel,
                                  conf.int    = TRUE,
                                  paired      = TRUE
                                )

          # Add row to main result table
          tableT$setRow( rowKey = key,
                         values = list(
                                        x        = xname,
                                        y        = yname,
                                        testName = "Wilcoxon W",
                                        t        = results$statistic,
                                        p        = results$p.value
                                      )
                       )
        } ## if (isParametric)

        # Test for normality; maybe add a warning marker to this pair of variables
        resultsS <- shapiro.test(x-y)
        if (isParametric) {
          if ( resultsS$p.value < 1-self$options$confLevel ) {
            tableT$addSymbol( rowKey=key, col=1, symbol='\u2021' )
            tableT$addSymbol( rowKey=key, col=2, symbol='\u2021' )
            warn_nn_ct <- warn_nn_ct + 1
          # } else {
          #   # If normality is not violated, possibly add a marker for significant result:
          #   if ( results$p.value < 1-self$options$confLevel ) {
          #     tableT$addSymbol( rowKey=key, col=1, symbol='*' )
          #      tableT$addSymbol( rowKey=key, col=2, symbol='*' )
          #      note_sig_ct <- note_sig_ct + 1
          #   }
          }
        }

        if ( self$options$normalityTest ) {
          # If user asked for display of normality test, add a line to normality test table:
          tableS$setRow( rowKey = key,
                         values = list(
                                        x = xname,
                                        y = yname,
                                        w = resultsS$statistic,
                                        p = resultsS$p.value
                                      )
                       )
        }
      }

      if ( warn_nn_ct > 0  ) tableT$setNote( 'normalWarning',
                                             paste( 'For pairs of variables marked \u2021 the assumption of normality of the differences may be violated.',
                                                    'Recommendation: use a nonparametric test (Wilcoxon Signed Rank Test).'
                                                  ),
                                             init=FALSE
                                           )
      # if ( note_sig_ct > 0 ) tableT$setNote( 'significantNote',
      #                                        paste( 'For pairs of variables marked * the calculated p-value is smaller than the chosen significance level.',
      #                                               'Therefore we reject the null hypothesis for each variable marked with *.'
      #                                             ),
      #                                        init=FALSE
      #                                      )

    } ## end .run

  )
)
