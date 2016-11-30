context("ColumnDef_constructor")

test_that("ColumnDef_object_is_created",{
  columnDef <- ColumnDef(columnName="myColumn",
                         categories=factor(c("Z","Y","X")),
                         type="categorical",
                         displayName="myDisplayName")
  
  expect_equal(class(columnDef)[1],"ColumnDef")
})

test_that("error_when_type_not_logical_categorical_or_numeric",{
  
  expect_error(ColumnDef(columnName="myColumn",
                         type="invalidType",
                         displayName="myDisplayName"))
  
  expect_error(ColumnDef(columnName="myColumn",
                         type=c("logical","categorical"),
                         displayName="myDisplayName"))
  
})

test_that("error_when_type_categorical_and_no_categories_given",{
  expect_error(ColumnDef(columnName="myColumn",
                         type="categorical",
                         displayName="myDisplayName",
                         unit="days"))
})

test_that("warning_when_include_categories_and_type_not_categorical",{
  
  expect_warning(ColumnDef(columnName="myColumn",
                         type="numeric",
                         displayName="myDisplayName",
                         categories=factor(c("c1","c2"))))
  
})

test_that("message_when_categories_not_factor",{
  
  expect_message(ColumnDef(columnName="myColumn",
                           type="categorical",
                           displayName="myDisplayName",
                           categories=c("c1","c2")))

})

test_that("warning_if_unit_not_character",{
  expect_warning(ColumnDef(columnName="myColumn",
            type="numeric",
            displayName="myDisplayName",
            unit=5))
})

test_that("error_if_columnName_length_not_1",{
  expect_error(ColumnDef(columnName=c("myColumn","c2"),
                           type="numeric",
                           displayName="myDisplayName"))
})

test_that("displayName_is_columnName_if_value_not_given",{
  columnDef <- ColumnDef(columnName="myColumn",
                         type="numeric")
  
  expect_equal(columnDef@displayName,"myColumn")
  
})

test_that("coercing_string_categories_to_factor_preserves_their_order",{
  
  expect_message(columnDef <- ColumnDef(columnName="myColumn",
                           categories=c("Z","Y","X"),
                           type="categorical",
                           displayName="myDisplayName"))
  
  expect_equal(levels(columnDef@categories),c("Z","Y","X"))
})

test_that("categories_given_as_factors_preserve_order",{
  columnDef <- ColumnDef(columnName="myColumn",
                         categories=factor(c("Z","Y","X")),
                         type="categorical",
                         displayName="myDisplayName")
  expect_equal(levels(columnDef@categories),levels(factor(c("Z","Y","X"))))
})

test_that("default_unit_is_empty_character",{
  columnDef <- ColumnDef(columnName="myColumn",
                         categories=factor(c("Z","Y","X")),
                         type="categorical",
                         displayName="myDisplayName")
  expect_equal(columnDef@unit,"")
})

test_that("non_categorical_columndefs_have_factor(NULL)_for_categories",{
  columnDef <- ColumnDef(columnName="myColumn",
                         type="logical")
  expect_equal(columnDef@categories,factor(NULL))
  
})