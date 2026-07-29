setwd(choose.dir()) #set the environment


library(tidyverse)
library(dplyr)
library(pracma)
library(ggforce)
library(purrr)

#### Load file (Spike List) - only set to do one file at a time ####

MEAdata <- read.csv(choose.files(default = "", caption = "Select files",
                                 multi = TRUE, filters = Filters,
                                 index = nrow(Filters)))

Filename <- MEAdata[29,2]           #Row location may vary with Axion updates
Filename <- gsub("\\(000\\)\\.raw", "", Filename)
Barcode <- MEAdata[19,2]            #Row location may vary with Axion updates

MEAdata <- MEAdata %>%
  select(3:5) %>%
  dplyr::mutate(Filename = Filename, Barcode = Barcode)

colnames(MEAdata)<- c("Time_s", "Electrode", "Amplitude_mV", "Filename", "Barcode")

tmax = 600 #time of recording in seconds

#Upload Platemap info

wellinfo<- read.csv(choose.files(default = "", caption = "Select files",
                                 multi = TRUE, filters = Filters,
                                 index = nrow(Filters)))#####add in well inf


##################Get just spikes data#########################

MEAdata$Amplitude_mV<-as.numeric(MEAdata$Amplitude_mV)
MEAdata$Time_s<-as.numeric(MEAdata$Time_s)


MEAdata <- MEAdata%>%
  filter(Amplitude_mV!= "NA") %>%
  separate(Electrode, into = c("Well", "Electrode"), sep= "_") 

MEAdata <- merge(MEAdata, wellinfo, by.x = "Well", by.y = "Well", all = FALSE, all.x = TRUE, all.y = TRUE)%>%
  filter(Active == "TRUE")

MEAdata<-arrange(MEAdata, Well, Time_s)


Spike_Stats <- MEAdata %>%
  dplyr::group_by(Filename,Treatment, Well, Electrode)%>%
  dplyr::mutate(TotalSpikesE = n()) %>%
  dplyr::ungroup()%>%
  dplyr::group_by(Filename,Well)%>%
  dplyr::mutate(TotalSpikesWell = n(), MedSpikesE = median(TotalSpikesE), 
                SpikeRateE = TotalSpikesE/tmax, MedSpikeRateE = median(SpikeRateE), 
                maxSpikerateE= max(SpikeRateE), SpikeRateWell = TotalSpikesWell/tmax, ActiveElec = n_distinct(Electrode), medAmplitude_mV = median (Amplitude_mV), meanAmplitude_mV = mean(Amplitude_mV))

Spike_Stats2<-Spike_Stats%>%
  select(-Time_s, -Electrode, -Amplitude_mV, -Well.Coloring) %>%
  distinct(Filename, Treatment, Well, .keep_all=TRUE)

csv_name <- paste(Barcode,"_", Filename, "_Spike_stats.csv", sep = "")

write.csv(Spike_Stats2, csv_name)

Total_Spikes <- Spike_Stats %>%
  select(Filename,Treatment, Well, TotalSpikesWell, ActiveElec) %>%
  distinct(Filename, Treatment, Well, .keep_all=TRUE)

##### Bin spikes, find max ASDR and SB threshold #####

# Step 1: Create a complete set of bins for each Filename and Well
complete_bins <- MEAdata %>%
  dplyr::group_by(Filename, Treatment, Well) %>%
  tidyr::expand(bin = cut(seq(0, 600, 0.2), seq(0, 600, 0.2))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(bin = as.character(bin))

# Step 2: Create your actual binned spike counts
Binned_Spikes <- MEAdata %>% 
  dplyr::group_by(Filename, Treatment, Well) %>%
  dplyr::mutate(bin = cut(Time_s, seq(0, 600, 0.2))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Filename, Treatment, Well, bin) %>%
  dplyr::summarize(count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(bin = as.character(bin))

# Step 3: Join complete bins with actual binned spikes and fill missing counts with zeros
Binned_Spikes <- complete_bins %>%
  dplyr::left_join(Binned_Spikes, by = c("Filename", "Well", "Treatment", "bin")) %>%
  dplyr::mutate(count = replace_na(count, 0))

# Step 4: Separate bin into start_time and end_time
Binned_Spikes <- Binned_Spikes %>%
  tidyr::separate(bin, c("start_time", "end_time"), sep = ",", remove = FALSE)%>%
  dplyr::mutate(end_time = gsub("]", "", end_time)) %>%
  dplyr::mutate(start_time = gsub("(", "", start_time, fixed=TRUE))

# Calculate Max ASDR and set 40% theshold
Well_MaxASDR <- Binned_Spikes %>%
  dplyr::group_by(Filename, Treatment, Well) %>%
  dplyr::mutate(maxASDR = as.numeric(max(count))) %>%
  dplyr::mutate(threshold = maxASDR*0.4)               ## Threshold set at 40% of MaxASDR

#### ASDR plots ####

Well_MaxASDR$end_time <- as.numeric(Well_MaxASDR$end_time)

split_ASDR <- split(Well_MaxASDR, Well_MaxASDR$Well)

ASDR_plot <- function(Well_MaxASDR) {
  p <- ggplot(Well_MaxASDR, aes(end_time, count)) +
    geom_line(stat = "identity", linejoin = "round", size = 0.1) +
    geom_hline(aes(yintercept = mean(threshold)), color = "red") +
    xlab("Time (s)")+
    ylab("Spike Count")+
    scale_x_continuous(limits = c(0, 600), breaks = seq(0, 600, 100))
  
  well_name <- unique(Well_MaxASDR$Well)
  Treatment <-unique(Well_MaxASDR$Treatment)
  ggsave(paste0(Barcode, "_",Filename, "_", well_name,"_",Treatment, "_ASDR_plot.tiff"), plot = p, height = 3, width = 10)
}

lapply(split_ASDR, ASDR_plot)


#### Calculate Synchronised Bust Interval, duration and number of spikes per SB ####

Sync_bursts <- Well_MaxASDR %>%                        ## Filter 200ms bins so that only the ones above 40% of maxASDR remain --> sync bursts ##
  filter(count>threshold
  )

Sync_bursts$end_time<-as.numeric(Sync_bursts$end_time)
Sync_bursts$start_time<-as.numeric(Sync_bursts$start_time)
Sync_bursts$count<-as.numeric(Sync_bursts$count)


SB_interval <- Sync_bursts %>%                     ## Calculate interval between sync bursts
  ungroup()%>%
  dplyr::group_by(Filename, Treatment, Well)%>%
  mutate(interval = start_time-lag(end_time))

Grouped_SB <- SB_interval %>%              ## will only group up to 5 bins (1 second)
  dplyr::group_by(Filename, Treatment, Well)%>%
  mutate(Spike_count = ifelse(interval < 0.4 & !is.na(interval), count + lag(count), count)) %>%       ##if interval is < 0.4 sec add spike count from previous row
  mutate(Spike_count = ifelse(interval < 0.4 & !is.na(interval), count + lag(Spike_count), Spike_count)) %>%
  mutate(Spike_count = ifelse(interval < 0.4 & !is.na(interval), count + lag(Spike_count), Spike_count)) %>%
  mutate(Spike_count = ifelse(interval < 0.4 & !is.na(interval), count + lag(Spike_count), Spike_count)) %>%  #repeat lag count so can get up to 5 bins 
  mutate(Spike_count = ifelse(interval < 0.4 & !is.na(interval), count + lag(Spike_count), Spike_count)) %>%  #combined together if interval is < 0.4
  mutate(start_time = ifelse(interval < 0.4 & !is.na(interval), lag(start_time), start_time)) %>%      ##if interval is < 0.4 sec change start time to match that of previous row
  mutate(start_time = ifelse(interval < 0.4 & !is.na(interval), lag(start_time), start_time)) %>%
  mutate(start_time = ifelse(interval < 0.4 & !is.na(interval), lag(start_time), start_time)) %>%
  mutate(start_time = ifelse(interval < 0.4 & !is.na(interval), lag(start_time), start_time)) %>%
  mutate(start_time = ifelse(interval < 0.4 & !is.na(interval), lag(start_time), start_time)) %>%     ## repeat lag 5 times so up to 4 bins can have same start time
  mutate(interval = ifelse(interval < 0.4, lag(interval), interval)) %>%                 ##if interval is < 0.4 sec change it so it matches previous row
  mutate(interval = ifelse(interval < 0.4, lag(interval), interval)) %>%
  mutate(interval = ifelse(interval < 0.4, lag(interval), interval)) %>%
  mutate(interval = ifelse(interval < 0.4, lag(interval), interval)) %>%
  mutate(interval = ifelse(interval < 0.4, lag(interval), interval))                    ## repeat lag 5 times

Grouped_SB2 <- Grouped_SB %>%
  dplyr::group_by(Filename, Treatment, Well)%>%
  filter(start_time != lead(start_time)) %>%   ## if start time matches next row, delete row
  drop_na(start_time) %>%
  mutate(SB_duration = end_time - start_time) %>%
  mutate(medianSBinterval = median(interval, na.rm = TRUE), meanSBinterval = mean(interval, na.rm = TRUE), 
         sd_SBinterval = sd(interval, na.rm = TRUE), SBinterval_CoV = sd_SBinterval/meanSBinterval) %>%          ## average SB interval
  mutate(medianSBduration = median(SB_duration, na.rm = TRUE), meanSBduration = mean(SB_duration, na.rm = TRUE),        ## average SB duration    
         no_SB = n(), median_spike_SB = median(Spike_count, na.rm = TRUE), mean_spike_SB = mean(Spike_count, na.rm = TRUE),         ##no of SBs, average no of spikes per SB
         max_spike_SB = max(Spike_count, na.rm = TRUE),                ## max no of spikes per SB - like max ASDR but combined from multiple 200ms bins
         total_spike_SB = sum(Spike_count, na.rm = TRUE))                ## total spikes in SBs - can be used to calculate % of spikes in SBs


## Simplify data frame to just wells, add total spikes data and calculate % spikes in bursts 
## and save output measures as csv

SB_stats <- Grouped_SB2 %>%
  select(-bin, -start_time, -end_time, -count, -interval, -Spike_count) %>%
  distinct(Filename, Treatment, Well, .keep_all=TRUE)

SB_stats <- SB_stats %>%
  left_join(Total_Spikes, by = join_by("Filename", "Treatment", "Well") ) %>%
  mutate(percent_spike_SB = total_spike_SB/TotalSpikesWell * 100)

csv_name <- paste(Barcode,"_", Filename, "_SB_stats.csv", sep = "")

write.csv(SB_stats, csv_name)

#### Spike Plots ####

MEAdata$Electrode <- factor(MEAdata$Electrode, levels = c("11", "12", "13", "14", "21", "22","23", "24", "31", "32", "33", "34", "41", "42", "43", "44"))

Split_Spike <- split(MEAdata, MEAdata$Well)

Spike_plot <- function(MEAdata) {
  p <- ggplot(MEAdata, aes(x = Time_s)) +
    geom_vline(aes(xintercept = Time_s)) +
    coord_cartesian(ylim = c(0,1))+
    facet_wrap(~Electrode, ncol = 1, strip.position = 'left', drop = FALSE) +
    scale_y_discrete(limits = c("11", "12", "13", "14", "21", "22","23", "24", "31", "32", "33", "34", "41", "42", "43", "44", "C", "B", "A")) +
    theme(axis.title.y       = element_blank(),
          axis.text.y        = element_blank(),
          axis.ticks.y       = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank())+
    xlab("Time (s)")+
    scale_x_continuous(limits = c(0, 600), breaks = seq(0, 600, 100))
  
  well_name <- unique(MEAdata$Well)
  Treatment <-unique(MEAdata$Treatment)
  
  ggsave(paste0(Barcode, "_", Filename, "_",  well_name,"_",Treatment, "_Spike_plot.tiff"), plot = p, height = 3, width = 10)
}

lapply(Split_Spike, Spike_plot)
