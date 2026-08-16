#THIS FUNCTION FILLS THE TAXONOMY OUTPUT FROM KRAKEN, WHICH HAS SOME NA SPACES, AS WELL AS "NA" STRINGS ("NA" BETWEEN TAXONOMIC RANKS)
fill_taxonomy <- function(taxonomy_data) {
  
  # Loop through each row of taxonomy_data
  for (i in 1:nrow(taxonomy_data)) {
    last_filled <- NULL  # # last_filled is a variable that stores the most recent non-NA value as we move through each row. It gets reset at the start of each row.
    
    # First loop: left to right to fill actual NAs
    for (j in 1:ncol(taxonomy_data)) {#This loop iterates through every column (j) in a given row (i). j starts at the first column\
      # (1 - leftmost column Domain) and increases (1 step at a time) until it reaches the last column (ncol(taxonomy_data) - Species).
      if (is.na(taxonomy_data[i, j])) { #This checks if the current cell (taxonomy_data[i, j]) contains an actual missing value (NA)
        # Handle actual NA (missing values): fill with 'unclassified' + last_filled
        #If last_filled is NULL (in case this is the first NA in the row), the value is filled with "unclassified". 
        #Otherwise if last_filled is not NULL (case where there is already a classification in the column that came before), 
        # the value is filled with "unclassified" and the most recent classification (stored in last_filled), forming a string like "unclassified [last classification]".
        taxonomy_data[i, j] <- ifelse(is.null(last_filled), "unclassified", paste("unclassified", last_filled)) 
      } else {
        last_filled <- taxonomy_data[i, j]  # If the current value is not NA (missing), it updates last_filled to this classification. This means that in the next iteration (as the loop moves right), this will be the latest known classification for filling unclassified fields.
      }
    }
    
    # Second loop: right to left to fill string "NA"s
    last_filled <- NULL  # Reset last_filled before the second loop
    for (j in ncol(taxonomy_data):1) {#This loop moves from the last column (Species) to the first (Kingdom) of every (i) given row. j starts at the rightmost column index (ncol(taxonomy_data) - Species) and decreases (1 step at a time) until it reaches the leftmost column (1 - Domain).
      if (taxonomy_data[i, j] == "NA") { ##Checks if the current cell contains "NA" (the string, not an actual missing value)
        # Look for the most recent valid value to the left. This nested for loop starts from the column immediately to the left of the 
        # current column (j-1) and moves towards the first column (1). It checks each value to see if it's not "NA". 
        #Once a valid (non-"NA") value is found, it is stored in correct_value, and the loop breaks (ends) because 
        #we've found the most recent valid value for this row.
        for (k in (j-1):1) {
          if (taxonomy_data[i, k] != "NA") {
            correct_value <- taxonomy_data[i, k]
            break
          }
        }
        
        # If a valid value is found (i.e., correct_value variable exists), replace "NA" with "unknown <correct_value> [<last_filled>]"
        if (exists("correct_value")) {
          taxonomy_data[i, j] <- paste0("unknown ", correct_value, " [", last_filled, "]")
        } else {
          taxonomy_data[i, j] <- paste0("unknown [", last_filled, "]") #If correct_value wasn't found (no valid non-"NA" value to the left), we just replace "NA" with "unknown [<last_filled>]"
        }
        
        # Clear correct_value after assignment
        rm(correct_value) 
        #After using correct_value, we remove it from memory with rm(correct_value) to avoid any issues when processing the next column. 
        #This makes sure the correct_value from the previous loop does not interfere with the next "NA" in the same row.
      } else if (!is.na(taxonomy_data[i, j])) { 
        last_filled <- taxonomy_data[i, j]  # If the value in the current column is not "NA" (it is a valid classification), 
        #we update last_filled with this valid value. This allows us to keep track of the most recent valid taxonomic classification,
        #which will be used in the next iterations when constructing the output for "NA" values.
      }
    }
  }
  
  return(taxonomy_data)
}