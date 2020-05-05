# Total USA: transform date and cases

covid <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv")
covid$date <- as.Date(covid$date)
mindate <- as.numeric(min(covid$date))
covid.total <- covid %>% group_by(date) %>% summarize(cases = sum(cases)) %>% 
  mutate(t = (as.numeric(date) - mindate)/100, y = cases/ max(cases)) 
covid <- covid.total



# Read data and transform date and cases
# change state name as desired

covid <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv")
covid$date <- as.Date(covid$date)
covid <- covid %>% filter(state == "South Carolina") 
mindate <- as.numeric(min(covid$date))
covid <- covid %>% mutate(t = (as.numeric(date) - mindate)/100, y = cases / max(cases))



# check if data was read and see plot

plot(covid$t, covid$y)
plot(covid$date, covid$cases)