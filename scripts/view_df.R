library(readr)

###A very powerful debugging function, opens up a current dataframe in visidata, inspired by a Reddit post lol.
view_df = function(df){
  df_string = format_tsv(df)
  system2("vd", input = df_string)
}
