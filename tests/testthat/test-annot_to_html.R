test_that("annot_to_html()", {
  expect_true(
    annot_to_html("[M+H]+_13C2", "C6H12O2", "M") == 
      "C<sub>4</sub><sup>13</sup>C<sub>2</sub>H<sub>13</sub>O<sub>2</sub><sup>+</sup>"
  )
  expect_true(
    annot_to_html("C6H13O2-H2O+_13C2_18O", "", "M") == 
      "C<sub>4</sub><sup>13</sup>C<sub>2</sub>H<sub>11</sub><sup>18</sup>O<sup>+</sup>"
  )
  expect_true(
    annot_to_html("C6H13O2-H2O\u2022_13C2_18O", "", "M") == 
      "C<sub>4</sub><sup>13</sup>C<sub>2</sub>H<sub>11</sub><sup>18</sup>O<sup>-</sup><sup>&#x2022</sup>"
  )
})

