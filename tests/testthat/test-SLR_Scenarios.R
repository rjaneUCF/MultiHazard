
test_that("SLR_Scenarios works", {


  result <-  SLR_Scenarios(SeaLevelRise=0.8, Scenario="Compact", Unit = "m", Year=2022,
                           Location="Key West")

  # Checking type of output
  expect_type(result, 'list')
  expect_named(result,c("High", "Intermediate", "Low"))

})


test_that("Invalid inputs", {
  expect_error(SLR_Scenarios(), "Error: SeaLevelRise missing")

  expect_error(SLR_Scenarios(SeaLevelRise = 0.5, Scenario = "BadScenario"),
               "Error: Invalid Scenario")

  expect_error(SLR_Scenarios(SeaLevelRise = 0.5, Unit = "feet"),
               "Error: Unit input must be m or Inches")

  expect_error(SLR_Scenarios(SeaLevelRise = 0.5, Unit = "m", Year= 1950),
               "Error: Invalid Year input")

  expect_error(SLR_Scenarios(SeaLevelRise = 0.5, Unit = "m", Location = 2030),
               "Error: Location name invalid")

})


test_that("Function is deterministic", {

 result1 <-  SLR_Scenarios(SeaLevelRise=0.8, Scenario="Compact", Unit = "m", Year=2022,
                           Location="Key West")

 result2 <-  SLR_Scenarios(SeaLevelRise=0.8, Scenario="Compact", Unit = "m", Year=2022,
                           Location="Key West")

expect_equal(result1, result2)

})

