theme1 = theme(
   # LABLES APPEARANCE
   panel.grid.major = element_blank(), 
   panel.grid.minor = element_blank(),
   panel.background = element_rect(fill = "transparent",colour = NA),
   plot.background = element_rect(fill = "transparent",colour = NA),
   plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
   axis.title.x = element_text(size=16, face="bold", colour = "black"),    
   axis.title.y = element_text(size=16, face="bold", colour = "black"),    
   axis.text.x = element_text(size=16, face="bold", colour = "black"), 
   axis.text.y = element_text(size=16, face="bold", colour = "black"),
   strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
   strip.text.y = element_text(size = 14, face="bold", colour = "black"),
   axis.line.x = element_line(color="black", size = 0.3),
   axis.line.y = element_line(color="black", size = 0.3),
   panel.border = element_rect(colour = "black", fill=NA, size=0.3)
   )


themegg = theme(
   # LABLES APPEARANCE
   panel.grid.major = element_blank(), 
   panel.grid.minor = element_blank(),
   panel.background = element_rect(fill = "transparent",colour = NA),
   plot.background = element_rect(fill = "transparent",colour = NA),
   plot.title = element_text(hjust = 0.5, size=14, face= "bold", colour= "black" ),
   axis.title.x = element_text(size=20, face="bold", colour = "black"),    
   axis.title.y = element_text(size=20, face="bold", colour = "black"),    
   axis.text.x = element_text(size=18, colour = "black"), 
   axis.text.y = element_text(size=18, colour = "black"),
   strip.text.x = element_text(size = 12, face="bold", colour = "black" ),
   strip.text.y = element_text(size = 12, face="bold", colour = "black"),
   strip.background =  element_rect(fill = "transparent",colour = NA),
   axis.line.x = element_line(color="black", linewidth = 0.2),
   axis.line.y = element_line(color="black", linewidth =  0.2),
   panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.2),
   legend.title=element_blank(),
   legend.text=element_text(size=12, face="bold", colour = "black"),
   legend.position="none"
)
