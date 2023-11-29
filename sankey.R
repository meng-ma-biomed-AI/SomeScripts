## libraries ----
library(ggalluvial)
library(ggrepel)
library(ggforce)
library(plotly)

## functions needed ----
cbindPad <- function(...){
  args <- list(...)
  n <- sapply(args,nrow)
  mx <- max(n)
  pad <- function(x, mx){
    if (nrow(x) < mx){
      nms <- colnames(x)
      padTemp <- matrix(NA, mx - nrow(x), ncol(x))
      colnames(padTemp) <- nms
      if (ncol(x)==0) {
        return(padTemp)
      } else {
        return(rbind(x,padTemp))
      }
    }
    else{
      return(x)
    }
  }
  rs <- lapply(args,pad,mx)
  return(do.call(cbind,rs))
}

## load data ----
sankey_data <- read.table('sankey_data.txt', sep = '\t', header = TRUE)

## data cleanup ----
# get unique IDs
pat_ids <- unique(sankey_data$person_id)
# blank df to hold data
combos_df <- data.frame(combo_ex = 'RVD')
# loop through to get data into necessary format
for (p in 1:length(pat_ids)) {
  combos_person <- sankey_data[sankey_data$person_id == pat_ids[p], ]
  person_unique_combos <- unique(combos_person$lot_drug)
  person_total <- data.frame(person_unique_combos)
  colnames(person_total) <- pat_ids[p]
  combos_df <- cbindPad(combos_df, person_total)
}
combos_df_t <- as.data.frame(t(combos_df))
combos_df_t <- combos_df_t[-1, ] # remove starter data

# format column names
colnames(combos_df_t) <- c('Treatment_1', 'Treatment_2', 'Treatment_3', 
                           'Treatment_4', 'Treatment_5', 'Treatment_6')

# subset to do Sankey on first 3 treatments
combos_df_t <- combos_df_t[, 1:3]

# assign drugs to categories
best_tx <- c('A', 'K', 'B', 'I')
better_tx <- c('J', 'C', 'F', 'E')

# loop through to re-assign
for (i in 1:nrow(combos_df_t)) {
  if (combos_df_t$Treatment_1[i] %in% best_tx) {
    combos_df_t$Treatment_1[i] <- 'best_tx'
  }
}

for (i in 1:nrow(combos_df_t)) {
  if (combos_df_t$Treatment_1[i] %in% better_tx) {
    combos_df_t$Treatment_1[i] <- 'better_tx'
  }
}

for (i in 1:nrow(combos_df_t)) {
  if (!is.na(combos_df_t$Treatment_1[i])) {
    if (combos_df_t$Treatment_1[i] != 'best_tx' & 
        combos_df_t$Treatment_1[i] != 'better_tx' & 
        combos_df_t$Treatment_1[i] != 'M') {
      combos_df_t$Treatment_1[i] <- 'good_tx'
    }
  }
}


for (i in 1:nrow(combos_df_t)) {
  if (combos_df_t$Treatment_2[i] %in% best_tx) {
    combos_df_t$Treatment_2[i] <- 'best_tx'
  }
}

for (i in 1:nrow(combos_df_t)) {
  if (combos_df_t$Treatment_2[i] %in% better_tx) {
    combos_df_t$Treatment_2[i] <- 'better_tx'
  }
}

for (i in 1:nrow(combos_df_t)) {
  if (!is.na(combos_df_t$Treatment_2[i])) {
    if (combos_df_t$Treatment_2[i] != 'best_tx' & 
        combos_df_t$Treatment_2[i] != 'better_tx' & 
        combos_df_t$Treatment_2[i] != 'M') {
      combos_df_t$Treatment_2[i] <- 'good_tx'
    }
  }
}

for (i in 1:nrow(combos_df_t)) {
  if (combos_df_t$Treatment_3[i] %in% best_tx) {
    combos_df_t$Treatment_3[i] <- 'best_tx'
  }
}

for (i in 1:nrow(combos_df_t)) {
  if (combos_df_t$Treatment_3[i] %in% better_tx) {
    combos_df_t$Treatment_3[i] <- 'better_tx'
  }
}

for (i in 1:nrow(combos_df_t)) {
  if (!is.na(combos_df_t$Treatment_3[i])) {
    if (combos_df_t$Treatment_3[i] != 'best_tx' & 
        combos_df_t$Treatment_3[i] != 'better_tx' & 
        combos_df_t$Treatment_3[i] != 'M') {
      combos_df_t$Treatment_3[i] <- 'good_tx'
    }
  }
}

for (i in 1:nrow(combos_df_t)) {
  for (j in 1:ncol(combos_df_t)) {
    if (is.na(combos_df_t[i,j])) {
      combos_df_t[i,j] <- 'no_tx'
    }
  }
}

## create and plot Sankeyy ----
df_sankey <- gather_set_data(combos_df_t, 1:3)
aw <- 0.03 # sets axis width
sp <- 0.2 # sets spacing
Unit <- rep(1, nrow(df_sankey)) # this is a dummy variable I couldn't plot without??

# plot it
p <- ggplot(df_sankey,
            aes(x = x, id = id, split = y, value = Unit, label = y)) +
  geom_parallel_sets(aes(fill = Treatment_1), alpha = 0.3, # fill needs to be adjusted for each level plotted
                     axis.width = aw, sep = sp, na.rm = FALSE) +
  geom_parallel_sets_axes(axis.width = aw, sep = sp, fill = c('black')) + # adjust manually
  theme_no_axes() + theme(legend.position = 'none') + 
  labs(x = '') + theme(axis.title.y=element_blank(),
                                                 axis.text.y=element_blank(),
                                                 axis.ticks.y=element_blank(), 
                                                 plot.title = element_text(face = 'bold'))


p + scale_fill_manual(values = c('dodgerblue', 'deeppink', 'seagreen3', 'darkorchid3')) + # add color after first level
  annotate(geom = 'text', x = 1, y = 1150, label = '1L', fontface = 'bold', size = 6) + 
  annotate(geom = 'text', x = 2, y = 1150, label = '2L', fontface = 'bold', size = 6) + 
  annotate(geom = 'text', x = 3, y = 1150, label = '3L', fontface = 'bold', size = 6)

p <- ggplot(df_sankey,
            aes(x = x, id = id, split = y, value = Unit, label = y)) +
  geom_parallel_sets(aes(fill = Treatment_2), alpha = 0.3, # fill needs to be adjusted for each level plotted
                     axis.width = aw, sep = sp, na.rm = FALSE) +
  geom_parallel_sets_axes(axis.width = aw, sep = sp, fill = c('black')) + # adjust manually
  theme_no_axes() + theme(legend.position = 'none') + 
  labs(x = '') + theme(axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(), 
                       plot.title = element_text(face = 'bold'))


p + scale_fill_manual(values = c('dodgerblue', 'deeppink', 'seagreen3', 'darkorchid3', 'gray')) + # add color after first level
  annotate(geom = 'text', x = 1, y = 1150, label = '1L', fontface = 'bold', size = 6) + 
  annotate(geom = 'text', x = 2, y = 1150, label = '2L', fontface = 'bold', size = 6) + 
  annotate(geom = 'text', x = 3, y = 1150, label = '3L', fontface = 'bold', size = 6)
