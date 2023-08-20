################################
# Implementation of Table Joins
################################

sort_merge_join <- function(x, table) {
  ox <- radixorderv(x)
  ot <- radixorderv(table)
  .Call(C_sort_merge_join, x, table, ox, ot)
}
