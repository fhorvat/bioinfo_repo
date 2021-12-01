library("rvest")
lego_movie <- read_html("http://www.imdb.com/title/tt1490017/")

lego_movie %>% 
  html_node("#title-overview-widget .small") %>%
  html_text() %>%
  as.numeric()

lego_movie %>%
  html_nodes(".itemprop .itemprop") %>%
  html_text()

lego_movie %>%
  html_nodes("table") %>%
  .[[2]] %>%
  html_table()

#tripadvisor demo
url <- "http://www.tripadvisor.com/Hotel_Review-g37209-d1762915-Reviews-JW_Marriott_Indianapolis-Indianapolis_Indiana.html"

reviews <- url %>%
  read_html() %>%
  html_nodes("#REVIEWS .innerBubble")

id <- reviews %>%
  html_node(".quote a") %>%
  html_attr("id")

quote <- reviews %>%
  html_node(".quote span") %>%
  html_text()

rating <- reviews %>%
  html_node(".rating .rating_s_fill") %>%
  html_attr("alt") %>%
  gsub(" of 5 stars", "", .) %>%
  as.integer()

date <- reviews %>%
  html_node(".rating .ratingDate") %>%
  html_attr("title") %>%
  strptime("%b %d, %Y") %>%
  as.POSIXct()

review <- reviews %>%
  html_node(".entry .partial_entry") %>%
  html_text()

data.frame(id, quote, rating, date, review, stringsAsFactors = FALSE) %>% View()