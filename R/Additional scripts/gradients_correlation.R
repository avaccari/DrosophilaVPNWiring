#just need to add r^2 values to the resulting plot (or any alternative measure)
#also is it possible to automate the titles for x and y axes? so that they would correspond to post1 and post2


con <- readRDS("data/hemibrain_con0.rds")

tmp0 <- con %>%
  filter(pre.type=="LC4", post.type %in% c("DNp02")) %>%
  group_by(pre.bodyID,post.type) %>%
  count() %>%
  ungroup()

tmp1 <- con %>%
  filter(pre.type=="LC4", post.type %in% c("DNp11")) %>%
  group_by(pre.bodyID, post.type) %>%
  count() %>%
  ungroup()

merge <- merge(tmp0,tmp1, by="pre.bodyID", all=TRUE)
merge[is.na(merge)] <- 0
merge<-merge%>%
  arrange(desc(n.x))


ggplot(merge, aes(x=n.x, y=n.y, col=pre.bodyID,))+
  geom_point(size = 3, col="steelblue")+
  geom_smooth(method="lm",formula= 'y ~ x', col="red", se=TRUE)+
  theme_classic()+
  ylab("LC4>DNp02 synapses")+
  xlab("LC4>DNp11 synapses")+
  theme(axis.text.x = element_text(face="bold", color="black", size=13, angle=0),
        axis.text.y = element_text(face="bold", color="black", size=15, angle=0))+
  theme(plot.title = element_text(face="bold", color="blue", size=15))+
  theme(axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15))
