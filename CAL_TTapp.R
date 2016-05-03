#################################################
# Title: Apps for Maximizing two response surfaces
# Author: Yongtao Cao
# Date: 05/03/2016
# Note: In order to run this code, you need to make
#       sure that all the following packages are 
#       installed AND updated to the latest version.
##################################################


library(shiny)
library(ggplot2)
library(gridExtra)
library(RColorBrewer) 
# Function to do Pareto comparison for maximizing 2 objectives

compare = function(newpt, curpf)
{
	g1 = round(newpt[1L], 4L) < round(curpf[,1L], 4L) 
	g2 = round(newpt[2L], 4L) < round(curpf[,2L], 4L) 
	
	ge1 = round(newpt[1L], 4L) <= round(curpf[,1L], 4L)
	ge2 = round(newpt[2L], 4L) <= round(curpf[,2L], 4L)
	
	l1 = round(newpt[1L], 4L) > round(curpf[,1L], 4L)
	l2 = round(newpt[2L], 4L) > round(curpf[,2L], 4L)
	
	le1 = round(newpt[1L], 4L) >= round(curpf[,1L], 4L)
	le2 = round(newpt[2L], 4L) >= round(curpf[,2L], 4L)
	
	eq1 = round(newpt[1L], 4L) == round(curpf[,1L], 4L)
	eq2 = round(newpt[2L], 4L) == round(curpf[,2L], 4L)
	
        
        cond1 = g1*ge2+g2*ge1 == 0
	cond2 = sum(l1*le2+l2*le1+eq1*eq2)
	
	newpf = curpf[cond1,]
	if(cond2 == 0)
	{
		newpf = rbind(newpf, newpt)
	}

        matrix(newpf, ncol = 2L)
}

# count function to see which point in the search space that is 
# on the Pareto front

cnt = function(A, B){
ind = rep(FALSE, nrow(A))
for (i in 1L:nrow(A)){
 for (j in 1L:nrow(B)){
  if(all(round(A[i,], 4L) == round(B[j,], 4L))){ 
        ind[i] = TRUE
} 
}
}
ind
}

sort.data.frame = function(tmpdata, tmpcol, or = 'FALSE')
{
	ndata = nrow(tmpdata)
	nsort = length(tmpcol)
	
	tmps = data.frame(array(0, c(ndata, 3L)))
        
	for (i in 1L:nsort)
	{
		tmps[,i] = tmpdata[c(tmpcol[i])]
	}
	if(or != "D")
	{
		index = order(tmps[,1L], tmps[,2L], tmps[,3L], decreasing = FALSE)
	}
	if(or == "D")
	{
		index = order(tmps[,1L], tmps[,2L], tmps[,3L], decreasing = TRUE)
	}
		
	tmpdata = tmpdata[index,]
	rownames(tmpdata) = c(1L:ndata)
	tmpdata
}

# Calculate 2D

vol2d = function(pfsc, r)
{
     pfsc[,1] = ifelse(pfsc[,1]==0, 1, (pfsc[,1]-min(pfsc[,1]))/(-min(pfsc[,1])))
     pfsc[,2] = ifelse(pfsc[,2]==0, 1, (pfsc[,2]-min(pfsc[,2]))/(-min(pfsc[,2])))

	if(nrow(pfsc) == 0L)
	{
		vol = 0L
	}else
	{
		temppf = sort.data.frame(pfsc, names(pfsc)[2L], or = "D")
		if(nrow(temppf) == 1L)
		{
			vol = sum((temppf[,1L] + r)*(temppf[,2L] + r))
		}else
		{
			vol = sum((temppf[,1L] + r)*(temppf[,2L] + r)) - sum((temppf[-(nrow(temppf)), 1L] + r)*(temppf[-1L, 2L] + r))	
		}
	}
	vol
}






canonrsm.v3 = function(statpt,ystat,eigratio,evecang,ngrid, target = 0) {          
stand.eigv1 = c(cos(evecang/180*pi),sin(evecang/180*pi))                     
stand.eigv2 = c(-sin(evecang/180*pi),cos(evecang/180*pi))                   
x1 = x2 = -ngrid:ngrid/ngrid
yval = matrix(999,ncol=2*ngrid+1,nrow=2*ngrid+1)
alldat = NULL
     for (i in 1:(2*ngrid+1)) {
      for (j in 1:(2*ngrid+1)) {
       if (target == 0){
       yval[i,j] = ystat + (-1)*eigratio*(t(stand.eigv1)%*%(c(x1[i],x2[j])-statpt))^2 +   
                     (-1)*(t(stand.eigv2)%*%(c(x1[i],x2[j])-statpt))^2}
       if (target > 0){ 
       yval[i,j] = abs(ystat + (-1)*eigratio*(t(stand.eigv1)%*%(c(x1[i],x2[j])-statpt))^2 +   
                     (-1)*(t(stand.eigv2)%*%(c(x1[i],x2[j])-statpt))^2 - target)}
       
       alldat = rbind(alldat,c(x1[i],x2[j],yval[i,j]))
}
}

a1 <- stand.eigv1[1]
b1 <- stand.eigv1[2]
c1 <- t(stand.eigv1)%*%statpt
a2 <- stand.eigv2[1]
b2 <- stand.eigv2[2]
c2 <- t(stand.eigv2)%*%statpt
beta0 <- ystat + (-1)*eigratio*c1^2 + (-1)*c2^2                             # for all betas
beta1 <- -2*(-1)*eigratio*a1*c1 -2*(-1)*a2*c2                               # replace eig1 with -1*eigratio
beta2 <- -2*(-1)*eigratio*b1*c1 -2*(-1)*b2*c2                               # replace eig2 with -1
beta11 <- -eigratio*a1^2 + (-1)*a2^2
beta22 <- -eigratio*b1^2 + (-1)*b2^2
beta12 <- -eigratio*2*a1*b1 + (-1)*2*a2*b2
return(list(c(beta0,beta1,beta2,beta11,beta22,beta12),alldat))
}


# Generate the search space

ui = (fluidPage(

   titlePanel("Multiple Optimal Design: Target VS Target"),

    sidebarLayout(
    sidebarPanel(
      sliderInput("M1X1", "X1 location for Y1:", min = -1, max = 1, value = 0.5, step = 0.2),
      sliderInput("M1X2", "X2 location for Y1:", min = -1, max = 1, value = -0.5, step = 0.2),
      sliderInput("Y1", "Maximum value of Y1:", min = 0, max = 250, value = 100, step = 1),
      sliderInput("YT1", "Target value of Y1:", min = 0, max = 250, value = 95, step = 1),
      sliderInput("er1", "Ratio of eigenvalues for Y1:", min = -10, max = 15, value = 10, step = 0.5),
      sliderInput("M1Ea", "Angle of max change for eigenvectors of Y1:", min = 0, max = 180, value = 45, step = 5),
      sliderInput("M2X1", "X1 location for Y2:", min = -1, max = 1, value = 0.5, step = 0.2),
      sliderInput("M2X2", "X2 location for Y2:", min = -1, max = 1, value = 0.5, step = 0.2),
      sliderInput("Y2", "Maximum value of Y2:", min = 0, max = 250, value = 50, step = 1),
      sliderInput("YT2", "Target value of Y2:", min = 0, max = 250, value = 48, step = 1),
      sliderInput("er2", "Ratio of eigenvalues for Y2:", min = -10, max = 15, value = 5, step = 0.5),
      sliderInput("M2Ea", "Angle of max change for eigenvectors of Y2:", min = 0, max = 180, value = 90, step = 5)
      
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("plot", height = 1000, width = 1000)
     
    )
  )
))

server = function(input, output) {

output$plot = renderPlot({
    
       temp.r1 = canonrsm.v3(c(input$M1X1,input$M1X2),input$Y1,input$er1,input$M1Ea,20, input$YT1) 
       temp.r2 = canonrsm.v3(c(input$M2X1,input$M2X2),input$Y2,input$er2,input$M2Ea,20, input$YT2)

       cdat1 = as.data.frame(temp.r1[[2]])
       cdat2 = as.data.frame(temp.r2[[2]])
       colnames(cdat1) = c("X1", "X2", "Dist1")
       colnames(cdat2) = c("X1", "X2", "Dist2")

    cont1 = ggplot(cdat1, aes(X1, X2, z = Dist1))+
            geom_tile(aes(fill = Dist1))+ 
            stat_contour(size = 1.5, color = "black")+
            scale_fill_gradientn(colours=brewer.pal(6,"YlOrRd"))+
            theme(axis.text = element_text(colour = 'black', size=16),
                  axis.title = element_text(size = 16), legend.position="bottom")
            
    cont2 = ggplot(cdat2, aes(X1, X2, z = Dist2))+
            geom_tile(aes(fill = Dist2))+ 
            stat_contour(size = 1.5, color = "black")+
            scale_fill_gradientn(colours=brewer.pal(6,"YlOrRd"))+
            theme(axis.text = element_text(colour = 'black', size=16),
                  axis.title = element_text(size = 16), legend.position="bottom")
       
       commat = cbind(temp.r1[[2]],temp.r2[[2]][,3])

       # Generate the Pareto fronts

       pf0 = cbind(temp.r1[[2]][,3],temp.r2[[2]][,3])

        newpf = pf0[1,,drop = FALSE]
	 for(i in 1L:nrow(pf0)){
           temppf = compare(matrix(pf0[i,], ncol = 2L), newpf)
            newpf = temppf
         }

# Get the the Position of PF in the search space

  ind = cnt(pf0, temppf)

# Get data to plot

     pdat = as.data.frame(commat[ind,,drop=FALSE])
     colnames(pdat) = c("X1", "X2", "Dist1", "Dist2")
     vdat = data.frame(-pdat[,3], -pdat[,4])
     volume = round(vol2d(vdat, r = 0),4)
     pdes = ggplot(pdat,aes(x=X1, y=X2))+
            theme(axis.text = element_text(colour = 'black', size=16),
                  axis.title = element_text(size = 16), legend.position="none")+
            geom_point(shape = 19, colour = 'red', size = (pdat$Dist1 - min(pdat$Dist1))*5/(max(pdat$Dist1)-min(pdat$Dist1)+1)+3)+
            scale_x_continuous(limits = c(-1, 1), breaks=seq(-1, 1, by = 0.2))+
            scale_y_continuous(limits = c(-1, 1), breaks=seq(-1, 1, by = 0.1))

      ppf = ggplot(pdat,aes(x=Dist1, y=Dist2))+
            theme(axis.text = element_text(colour = 'black', size=16),
                  axis.title = element_text(size = 16), legend.position="none")+
            geom_step(color='blue', direction = "hv", size = 1)+
            geom_point(colour = 'red', size = (pdat$Dist1 - min(pdat$Dist1))*5/(max(pdat$Dist1)-min(pdat$Dist1)+1)+3)+
            geom_point(x = 0, y = 0, shape = 42, size = 16)+
            scale_x_continuous(limits = c(0, max(pdat$Dist1)))+
            scale_y_continuous(limits = c(0, max(pdat$Dist2)))+
	    geom_text(aes(x = max(pdat$Dist1)-diff(range(pdat$Dist1))/4, 
                          y = max(pdat$Dist2)-diff(range(pdat$Dist2))/6, 
                          label=paste("Hypervolume: ", volume)), size = 8)

   grid.arrange(cont1, cont2, pdes, ppf, ncol = 2)
   })

}

shinyApp(ui = ui, server = server)



