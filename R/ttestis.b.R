# Implementation of the independent two sample test
#
# This implementation is NOT optimized.
# * It could make use of jamovi's clearWith feature (but it is a bit tricky to get this right for all cases).
# * Shapiro-Wilks and Levene tests could be cached (per variable), because they are independent of all other settings.
# * t-test/Mann-Whitney test results could be cached (per variable, confLevel, mean, and alternative).
# But all the tests are fast, so optimisation probably does not pay off.
#

ttestISClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
  "ttestISClass",
  inherit = ttestISBase,
  private = list(

    .init = function( ) {
      # Do some initialisations on the result tables that could not be done statically via ttestis.a.yaml

      isParametric   <- self$options$testType == 'parametric'

      tableT         <- if (isParametric) self$results$ttestParIS else self$results$ttestNonparIS
      tableL         <- self$results$homogeneityIS
      tableS         <- self$results$normalityIS

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

        tableT$setNote( 'nullhyp', 'Null hypothesis: the difference of means between the groups in the population is equal to 0', init=TRUE )

        descr <- paste( 'Alternative hypothesis: the difference of means between the groups in the population is',
                        if ( self$options$alternative=='two.sided' ) 'not equal to 0' else paste( self$options$alternative, 'than 0' )
                      )
        tableT$setNote( 'althyp', descr, init=TRUE )
      }

      descr <- paste0( 'Null hypothesis (H0): the distributions in the groups have equal variances.\n',
                       'If the p-value is smaller than the significance level, we reject this H0.'
                     )
      tableL$setNote( 'homogeneity', descr, init=TRUE )

      descr <- paste0( 'Null hypothesis (H0): the data for each group in the population come from a normal (i.e. Gaussian) distribution.\n',
                       'If the p-value is smaller than the significance level, we reject this H0.'
                     )
      tableS$setNote( 'normality', descr, init=TRUE )

      self$results$ttestParIS$setVisible(isParametric)
      self$results$ttestNonparIS$setVisible(!isParametric)
      tableL$setVisible(isParametric)
      tableS$setVisible( self$options$normalityTest )

    }, ## end .init

    .run = function( ) {

      # `self$data`    contains the data
      # `self$options` contains the options
      # `self$results` contains the results object (to be populated here)

      data   <- self$data

      isParametric <- self$options$testType == 'parametric'
      tableT       <- if (isParametric) self$results$ttestParIS else self$results$ttestNonparIS
      tableL       <- self$results$homogeneityIS
      tableS       <- self$results$normalityIS

      grp          <- self$options$group
      if ( is.null(grp) ) return( )

      grplevels    <- levels( data[,grp] )
      if ( length(grplevels) != 2 ) {
        msg = paste( 'Grouping variable', grp, 'must have exactly 2 different values.',
                     'Consider going to Data -> Filter to include only 2 groups.'
                   )
        jmvcore::reject(msg)
      }

      rowTNo       <- 0
      rowSNo       <- 0
      warn_nn_ct   <- 0
      warn_var_ct  <- 0

      for ( dep in self$options$deps ) {

        # for all dependent variables, calculate the appropriate test, including Shapiro-Wilks, and add a line to the result table

        x <- jmvcore::toNumeric( data[ data[,grp]==grplevels[1], dep ] )
        y <- jmvcore::toNumeric( data[ data[,grp]==grplevels[2], dep ] )
        if ( is.factor(x) ) jmvcore::reject( paste( 'Cannot run test on grouping variable', dep ) )

        if (isParametric) {

          # Parametric case: both Student's t and Welch's t
          results  <- t.test( x, y,
                              alternative = self$options$alternative,
                              mu          = 0,
                              var.equal   = TRUE,
                              conf.level  = self$options$confLevel,
                              paired      = FALSE
                            )
          resultsW  <- t.test( x, y,
                               alternative = self$options$alternative,
                               mu          = 0,
                               var.equal   = FALSE,
                               conf.level  = self$options$confLevel,
                               paired      = FALSE
                             )

          df        <- results$parameter
          sided     <- if ( results$alternative=='two.sided' ) 2 else 1
          alpha     <- 1 - attr( results$conf.int, 'conf.level' )
          critval   <- qt( 1 - alpha/sided, df )

          dfW       <- resultsW$parameter
          alpha     <- 1 - attr( resultsW$conf.int, 'conf.level' )
          critvalW  <- qt( 1 - alpha/sided, dfW )

          if ( results$alternative=='less' ) {
            critval  <- -critval
            critvalW <- -critvalW
          }

          # Add results for Student's t-test to main table
          rowTNo     <- rowTNo + 1
          saveRowTNo <- rowTNo
          tableT$addRow( rowKey = rowTNo,
                         values = list(
                                    var      = dep,
                                    testName = "Student's t",
                                    t        = results$statistic,
                                    crit     = critval,
                                    df       = df,
                                    p        = results$p.value,
                                    mean1    = results$estimate[1],
                                    mean2    = results$estimate[2],
                                    meandiff = results$estimate[1] - results$estimate[2],
                                    cilower  = results$conf.int[1],
                                    ciupper  = results$conf.int[2]
                                  )
                       )

          # Add Welch results to main table:
          rowTNo <- rowTNo + 1
          tableT$addRow( rowKey = rowTNo,
                         values = list(
                                    var      = dep,
                                    testName = "Welch's t",
                                    t        = resultsW$statistic,
                                    crit     = critvalW,
                                    df       = dfW,
                                    p        = resultsW$p.value,
                                    mean1    = resultsW$estimate[1],
                                    mean2    = resultsW$estimate[2],
                                    meandiff = resultsW$estimate[1]-resultsW$estimate[2],
                                    cilower  = resultsW$conf.int[1],
                                    ciupper  = resultsW$conf.int[2]
                                  )
                       )

          # Calculate and add Levene test results (homogeneity of variances) to separate table
          resultsL = car::leveneTest( data[,dep], data[,grp], center='mean' )

          tableL$setRow( rowKey = dep,
                         values = list(
                                    var      = dep,
                                    f        = resultsL$'F value'[1],
                                    df       = resultsL$'Df'[1],
                                    df2      = resultsL$'Df'[2],
                                    p        = resultsL$'Pr(>F)'[1]
                                  )
                       )
          if ( resultsL$'Pr(>F)'[1] < 1-self$options$confLevel ) {
            tableT$addSymbol( rowNo=saveRowTNo, col=1, symbol='\u2020' )
            warn_var_ct <- warn_var_ct + 1
          }

        } else {

          # Non-parametric case: Mann-Whitney U test
          results <- wilcox.test(
                                   x, y,
                                   alternative = self$options$alternative,
                                   mu          = 0,
                                   digits.rank = 7,
                                   conf.level  = self$options$confLevel,
                                   conf.int    = TRUE,
                                   paired      = FALSE
                                )

          # Add results for Mann-Whitney to main table
          rowTNo <- rowTNo + 1
          saveRowTNo <- rowTNo
          tableT$addRow( rowKey = rowTNo,
                         values = list(
                                    var        = dep,
                                    testName   = 'Mann-Whitney U',
                                    t          = results$statistic,
                                    p          = results$p.value
                                  )
                       )

        }


        # Test for normality of either group; maybe add a warning marker to the variable
        resultsSx <- shapiro.test(x)
        resultsSy <- shapiro.test(y)

        if (isParametric) {
          if ( ( resultsSx$p.value < 1-self$options$confLevel ) || ( resultsSy$p.value < 1-self$options$confLevel ) ) {
            tableT$addSymbol( rowNo=rowTNo-1, col=1, symbol='\u2021' )
            warn_nn_ct <- warn_nn_ct + 1
          }
        }

        # Maybe display Shapiro-Wilk test for normality:
        if ( self$options$normalityTest ) {
          # If user asked for display of normality test, add a line to normality test table:
          rowSNo <- rowSNo + 1
          tableS$addRow( rowKey = rowSNo,
                         values = list(
                                        var   = dep,
                                        group = grplevels[1],
                                        w     = resultsSx$statistic,
                                        p     = resultsSx$p.value
                                      )
                       )
          rowSNo <- rowSNo + 1
          tableS$addRow( rowKey = rowSNo,
                         values = list(
                                        var   = dep,
                                        group = grplevels[2],
                                        w     = resultsSy$statistic,
                                        p     = resultsSy$p.value
                                      )
                       )
        }

      }

      if ( warn_nn_ct > 0  ) tableT$setNote( 'normalWarning',
                                             'For variables marked \u2021 the assumption of normality may be violated for at least one of the groups.',
                                             init=FALSE
                                           )
      if ( warn_var_ct > 0 ) tableT$setNote( 'varianceWarning',
                                             'For variables marked \u2020 the assumption of equal variances for the groups may be violated.',
                                             init=FALSE
                                           )

    } ## end .run

  )
)
