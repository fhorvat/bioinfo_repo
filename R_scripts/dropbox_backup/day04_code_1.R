setwd("C:/Users/Filip/Dropbox/Praksa bioinfo/Projekti/test/adventOfCode")
input1 <- read.table("day04_input_1.txt", sep = "", stringsAsFactors = F)

checkChecksum <- function(input_sub){
  
  # split string, get sector ID and checksum
  input_sub_split <- unlist(strsplit(input_sub, split = "-|\\[|\\]"))
  input_sub_ID <- input_sub_split[grep("[[:digit:]]", input_sub_split)]
  input_sub_checksum <- input_sub_split[length(input_sub_split)]
  
  # remove sector ID and checksum from string
  input_sub_split <- input_sub_split[-grep(input_sub_ID, input_sub_split)]
  input_sub_split <- input_sub_split[-grep(input_sub_checksum, input_sub_split)]
  
  # split to individual characters, get counts of characters, order by count
  input_sub_split <- unlist(strsplit(input_sub_split, split = ""))
  input_sub_split_table <- as.data.frame(table(input_sub_split), stringsAsFactors = F)
  input_sub_split_table <- input_sub_split_table[order(input_sub_split_table$Freq, decreasing = T), ]
  
  # get real checksum, compare with given checksum, if same return input sector ID
  input_sub_split_table_checksum <- paste0(input_sub_split_table$"input_sub_split"[1:5], collapse = "")
  if (input_sub_checksum == input_sub_split_table_checksum){
    return(as.integer(input_sub_ID))
  } else{
    return(0)
  }
}

all_sector_ID <- apply(input1, 1, checkChecksum)
sum(all_sectorID)


# part 2
letters_double <- c(letters, letters)

decryptRoomName <- function(input_sub){
  # split string, get sector ID and checksum
  input_sub_split_index <- grep("-", unlist(strsplit(input_sub, split = "")))
  input_sub_split <- unlist(strsplit(input_sub, split = "-|\\[|\\]"))
  input_sub_ID <- input_sub_split[grep("[[:digit:]]", input_sub_split)]
  input_sub_checksum <- input_sub_split[length(input_sub_split)]
  
  # remove sector ID and checksum from string
  input_sub_split <- input_sub_split[-grep(input_sub_ID, input_sub_split)]
  input_sub_split <- input_sub_split[-grep(input_sub_checksum, input_sub_split)]
  input_sub_split <- paste0(input_sub_split, collapse = "-")
  
  # split to individual characters
  input_sub_split <- unlist(strsplit(input_sub_split, split = ""))
  
  # decrypt = divide shift number with number of letters in english alphabet (=26), 
  # take reminder and shift letters by that reminder in vector which contains alphabet 
  # twice in a row
  input_sub_decrypt <- letters_double[match(input_sub_split, letters) + (as.integer(input_sub_ID) %% 26)]
  input_sub_decrypt[is.na(input_sub_decrypt)] <- " "
  input_sub_decrypt <- paste0(input_sub_decrypt, collapse = "")
  return(input_sub_decrypt)
}

all_decrypted_names <- apply(input1, 1, decryptRoomName)
all_decrypted_rooms_with_ID <- data.frame(room_name = all_decrypted_names, ID = all_sector_ID, stringsAsFactors = F) 
all_decrypted_rooms_with_ID <- all_decrypted_rooms_with_ID[all_decrypted_rooms_with_ID$ID != 0, ]
all_decrypted_rooms_with_ID[grep("northpole", all_decrypted_rooms_with_ID$room_name), "ID"]
