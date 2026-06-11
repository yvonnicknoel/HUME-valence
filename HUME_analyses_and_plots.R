emo = read.csv("emotions-2024-2026.csv")

# Gender + Age + Induction + 40 pre-induction responses + 40 post-induction responses
dim(emo)
# [1] 326  83

# Subjects' characteristics
subjects = emo[,1:3]

table(subjects$Induction)
# Angry    Fear   Happy Neutral     Sad 
#    62      72      67      61      64 

table(subjects$Gender)
#   F   M 
# 270  56 

summary(subjects$Age)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 18.00   19.00   20.00   22.92   23.00   59.00 

# Extract subtables: Emotion assessments before and after the induction
Pre = emo[,4:43]    # N x 40
Post = emo[,44:83]  # N x 40

# Homogeneize item labels (initially numbered)
emo.names = colnames(Pre)
emo.names = gsub("1$","",emo.names)
colnames(Pre) = emo.names
colnames(Post) = emo.names

# Pseudo bayesian correction for extreme responses on 0-100
nPre = (Pre+1) / 102
nPost = (Post+1) / 102

#--------------------- Figure 5: Illustration of the parameters within the Beta Unfolding (Noel, 2014) ------------------------------

# Load the BUM parameter estimation routine
source("BUM.R")

x11(height=8, width=13)
M.sim = BUM$new(2000, 6, delta=c(-2,0,0,2,0,0), lambda=c(2,1,0,2,3,3), tau=c(1,1,-1,1,1,2))
par(mfrow=c(2,3))
M.sim$Plot(type="true.model",items=1,mar=c(5,5,4,2),plot.data=FALSE,plot.locations=TRUE, main="Left locations")
M.sim$Plot(type="true.model",items=2,mar=c(5,5,4,2),plot.data=FALSE,plot.locations=TRUE, main="Low acceptability")
M.sim$Plot(type="true.model",items=3,mar=c(5,5,4,2),plot.data=FALSE,plot.locations=TRUE, plot.ci=TRUE, main="High dispersion")
M.sim$Plot(type="true.model",items=4,mar=c(5,5,4,2),plot.data=FALSE,plot.locations=TRUE, main="Right location")
M.sim$Plot(type="true.model",items=5,mar=c(5,5,4,2),plot.data=FALSE,plot.locations=TRUE, main="High acceptability")
M.sim$Plot(type="true.model",items=6,mar=c(5,5,4,2),plot.data=FALSE,plot.locations=TRUE, plot.ci=TRUE, main="Low dispersion")

dev.copy2pdf(file="Images/Fig.5-BUM-params.pdf")


#--------------------------------------- Beta Unfolding Analysis of Pre-induction data -------------------------------------

# We keep all items by default
item.sel = emo.names

Pre.sel = nPre[,item.sel]
Post.sel = nPost[,item.sel]

# Define those that are theoretically expected inverted
sel.inv = c("Passive","Bored","Asleep","Calm")

# Open graphical window: The iteratively adjusted ICC will be plotted during estimation
x11(width=19,height=10)
mfrow = c(5,8)

# Create the analysis object (load the dataset, specify reversed items)
M.pre = BUM$new(data=Pre.sel, inverted=sel.inv)

# Get MAP estimates of parameters
M.pre$Estimate(mfrow=mfrow)

# Cosmetic: Put positive emotion on the right
M.pre$theta = -M.pre$theta
M.pre$delta = -M.pre$delta

# Display ICCs with a bit more detail
M.pre$Plot(mfrow=mfrow, plot.ci=TRUE, plot.locations=TRUE, plot.smooth=TRUE)

dev.copy2pdf(file="Images/emotion40-BUM-pre.pdf")

# Parameter estimates
M.pre$Summary()

#                     DELTA      A.S.E       LAMBDA      A.S.E         TAU      A.S.E    OUTFIT     INFIT Misfit
# Depressed    -3.732286096 5.66966955  0.061479477 0.06958864 -1.76341978 5.63861629 1.0996439 1.0805877       
# Unhappy      -3.694759411 5.33334554  0.007946557 0.06918189 -1.53259671 5.30240766 1.0720214 1.2756688       
# Insecure     -3.652639947 5.08058483 -0.141146807 0.06846683 -2.01439926 5.04659871 0.9750269 1.1545801       
# Scared       -3.505742727 3.83241310 -0.122319140 0.06841052 -1.51078668 3.79972014 0.9508756 1.2932825       
# Disheartened -3.456304757 3.50916011 -0.166144009 0.06825466 -1.74203089 3.47522115 1.2664606 1.0997277       
# Panicked     -3.322563013 2.70135821 -0.103687228 0.06846141 -1.11132260 2.66961758 1.1691693 1.4602347    ***
# Frustrated   -2.533309025 0.68448534 -0.152841556 0.06837740 -0.93416463 0.65147354 1.1770950 1.0766758       
# Anxious      -2.490494015 0.61752468  0.023214079 0.07035806 -1.55193531 0.58405136 0.9505240 0.9370121       
# Angry        -2.440740538 0.56244481  0.032687740 0.06922094  0.21276031 0.53549790 1.2127396 1.7091602    ***
# Upset        -2.391814569 0.56664490 -0.244434537 0.06772774 -0.73029291 0.53299722 1.2514862 1.1220342       
# Aggressive   -2.300158046 0.44533561  0.090877605 0.06954762  0.63958294 0.42060315 1.1465708 1.5957610    ***
# Helpless     -2.238532949 0.45643745 -0.307935730 0.06762809 -1.03106080 0.42083448 1.5235727 1.0008631    ***
# Tense        -2.167222467 0.38002254  0.059633248 0.06988278 -0.76948965 0.34998672 0.9122714 0.9673702       
# Overwhelmed  -2.113136068 0.37936857 -0.273617455 0.06825614 -1.21997845 0.34350396 1.0694575 1.0129559       
# Nervous      -2.064022229 0.32310234  0.178353296 0.07076054 -0.78951492 0.29438534 0.8837993 0.8672249       
# Sad          -1.897872460 0.26983989  0.043289684 0.06956845 -0.30535001 0.24234820 1.0524494 0.9950739       
# Vulnerable   -1.538570185 0.19718855 -0.205240886 0.06853803 -0.69193307 0.16513905 1.1393986 0.9190271       
# Stressed     -1.528924073 0.17720119  0.201047476 0.07172150 -0.85755745 0.15014023 0.9786836 0.8484168       
# Calm         -0.604944284 0.10362141 -0.145476431 0.07041800 -0.90071293 0.08057910 1.0011483 0.9245251       
# Vigilant     -0.306703311 0.09433446 -0.129059181 0.07016509 -0.76913894 0.07647859 1.0710971 0.9936940       
# Elated       -0.071409154 0.10063999 -0.359406821 0.06711051 -0.02348304 0.08444371 1.0986667 1.0401379       
# Surprised    -0.005666605 0.10708504 -0.505115014 0.06603963  0.15413480 0.08822634 1.2771336 1.0920102       
# Unconcerned   0.023726621 0.10814574 -0.630237929 0.06584572 -0.53203804 0.08408642 1.3050433 0.9850824       
# InControl     0.447143852 0.08528350  0.203241261 0.07354257 -1.00128359 0.07570135 0.8936925 0.8646839       
# Empowered     0.447586654 0.09120277 -0.069283037 0.06993828 -0.61063503 0.08082998 0.9545267 0.9057702       
# Active        0.462708196 0.08924798  0.013667344 0.07104362 -0.77486649 0.07899516 0.9480477 0.9149149       
# Asleep        0.654642455 0.09027610  0.339974248 0.07991821 -1.65748539 0.07478547 1.0658784 1.0639988       
# Enthusiastic  0.673821627 0.08860292  0.282042686 0.07346268 -0.91804238 0.08010293 0.9005531 0.9007686       
# Energetic     0.675651851 0.08033185  0.614040855 0.07421586 -0.54043210 0.08029444 0.8771943 0.8677757       
# Passive       0.685406135 0.08842528  0.724652322 0.08405937 -1.70580192 0.07434218 1.0584438 1.0755492       
# Satisfied     0.729477595 0.08906910  0.308358260 0.07290015 -0.74983779 0.08230326 0.8880612 0.8805603       
# Happy         0.813301944 0.08450557  0.735311554 0.07564342 -0.76917120 0.08177393 0.8403972 0.8146576       
# Bored         0.835702312 0.09127279  1.254039270 0.08732247 -1.72292004 0.07781876 1.2398902 1.2798592       
# Content       0.839102907 0.06859527  1.609406056 0.07766469 -0.29080369 0.08056591 0.9417408 0.9338732       
# Joyful        0.865091090 0.09171218  0.477082399 0.07457844 -0.91273348 0.08405465 0.8221791 0.7896657       
# SelfAssured   1.068097019 0.10016974  0.500961037 0.07396675 -0.79111278 0.09201057 0.9080978 0.8454755       
# Confident     1.106789296 0.10282605  0.480451212 0.07379131 -0.79239941 0.09397643 0.9209252 0.8885404       
# Relaxed       1.349512383 0.11978783  0.576947952 0.07440611 -0.89573184 0.10764052 0.9193809 0.9278040       
# AtEase        1.372596864 0.12549857  0.469811651 0.07429282 -1.05930789 0.11045691 0.9319469 0.9177059       
# Serene        1.483485545 0.13229330  0.711271408 0.07543975 -1.02609164 0.11799042 0.8932876 0.8650684       

# NOTE: standard errors are huge on the negative valence side but this is easy to understand. Most subjects are
# on the positive side (see the joint person-item plot below), so the distance to the most negative emotions is large,
# hence the uncertainty. The emotion peak locations are well outside the observed range of person locations and the estimation
# routine tries its best to discriminate between slight slope variations to estimate distant locations.

#------------------------------- Figure 6: Selection of 15 illustrative ICC ------------------------------------
x11(height=10,width=7.5)

M.pre$Plot(mfrow=c(5,3), items=c("Depressed","Angry","Scared","Anxious","Stressed","Bored","InControl","Active","Calm","Enthusiastic","Joyful","Passive","Satisfied","Happy","Serene"),plot.ci=TRUE, plot.locations=TRUE, plot.smooth=TRUE)

dev.copy2pdf(file="Images/Fig.6-emotion15-BUM-pre.pdf")


#-------------------- Figure 7: Plot of emotion peak locations along the valence dimension, by domain ----------------------

# Pre-induction parameters are taken as reference parameters
params = M.pre$Summary(print.it=FALSE)$items

# Export item parameters and statistics to LaTeX
cat(M.pre$toLatex(params),"\n")

deltas = params[,1]
names(deltas) = rownames(params)
print(deltas)

#    Depressed      Unhappy     Insecure       Scared Disheartened     Panicked   Frustrated      Anxious        Angry        Upset 
# -3.732286096 -3.694759411 -3.652639947 -3.505742727 -3.456304757 -3.322563013 -2.533309025 -2.490494015 -2.440740538 -2.391814569 
#   Aggressive     Helpless        Tense  Overwhelmed      Nervous          Sad   Vulnerable     Stressed         Calm     Vigilant 
# -2.300158046 -2.238532949 -2.167222467 -2.113136068 -2.064022229 -1.897872460 -1.538570185 -1.528924073 -0.604944284 -0.306703311 
#       Elated    Surprised  Unconcerned    InControl    Empowered       Active       Asleep Enthusiastic    Energetic      Passive 
# -0.071409154 -0.005666605  0.023726621  0.447143852  0.447586654  0.462708196  0.654642455  0.673821627  0.675651851  0.685406135 
#    Satisfied        Happy        Bored      Content       Joyful  SelfAssured    Confident      Relaxed       AtEase       Serene 
#  0.729477595  0.813301944  0.835702312  0.839102907  0.865091090  1.068097019  1.106789296  1.349512383  1.372596864  1.483485545 

domains = list(
    Depression = c("Depressed", "Unhappy", "Disheartened", "Sad"),
    Fear = c("Insecure", "Scared", "Panicked"),
    Angriness = c("Angry", "Frustrated", "Aggressive", "Upset"),
    Vulnerability = c("Helpless", "Overwhelmed", "Vulnerable"),
    Stress = c("Stressed", "Anxious", "Tense", "Nervous", "Calm(-)"),
    "Arousal-Activity" = c("Active", "Surprised", "Vigilant", "Energetic", "Passive(-)", "Asleep(-)", "Unconcerned"),
    "Dominance-Control" = c("InControl", "Empowered", "SelfAssured", "Confident"),
    "Pleasure-Excitement" = c("Happy", "Joyful", "Enthusiastic", "Satisfied", "Content", "Elated", "Bored(-)"),
    Serenity =c("Serene", "Relaxed", "AtEase")
)

library(wordcloud)
x11(height=10, width=10)

# Note: Reversed items (marked as (-)) won't print (as strictly speaking, they have two locations)
nd = length(domains)
par(mar=c(5,2,2,2))
plot(NULL, type="n", xlab="Peak locations", ylab="", xlim=c(-6.5, 2.3), ylim=c(0, nd+1), main="", yaxt="n")
for(i in 1:nd) {
  text(-4.8,i,names(domains)[i], cex=1.1, col=gray(.35), adj=1)
  segments(-4.5,i,2,i)
  labels = domains[[i]]
  d = deltas[labels]
  k = order(d)
  d = d[k]
  labels = labels[k]
  for(j in 1:length(labels)) {
    text(d[j],i,"|")
    text(d[j],i+0.2+((0.2*(j-1)) %% 0.6),labels[j],cex=.8)
  }
}
rug(M.pre$theta, side=1)
text(0,0.3,"Pre-induction participants' distribution",col=gray(.35))

dev.copy2pdf(file="Images/Fig.7-item-locations.pdf")

# Note: Superimposed labels have been manually adjusted

#--------------------------------------- Analysis of post-induction data --------------------------------

# Group selectors
gAngry = subjects$Induction == "Angry"
gHappy = subjects$Induction == "Happy"
gNeutral = subjects$Induction == "Neutral"
gSad = subjects$Induction == "Sad"
gFear = subjects$Induction == "Fear"

# BUM analyses on separate groups, unconstrained
# The goal is to assess stability of item peak locations accross induction groups

M.post.neutral = BUM$new(data=Post.sel[gNeutral,], inverted=sel.inv)
M.post.neutral$Estimate(display=FALSE)

M.post.angry = BUM$new(data=Post.sel[gAngry,], inverted=sel.inv)
M.post.angry$Estimate(display=FALSE)

M.post.happy = BUM$new(data=Post.sel[gHappy,], inverted=sel.inv)
M.post.happy$Estimate(display=FALSE)

M.post.sad = BUM$new(data=Post.sel[gSad,], inverted=sel.inv)
M.post.sad$Estimate(display=FALSE)

M.post.fear = BUM$new(data=Post.sel[gFear,], inverted=sel.inv)
M.post.fear$Estimate(display=FALSE)

#--- Correlations between separate estimations of locations across measurement times and induction conditions
library(Hmisc)

locations = cbind(M.pre$delta,M.post.happy$delta,M.post.neutral$delta,M.post.sad$delta,M.post.angry$delta,M.post.fear$delta)
colnames(locations) = c("Pre","PostHappy","PostNeutral","PostSad","PostAngry","PostFear")
rcorr(locations)

#               Pre PostHappy PostNeutral PostSad PostAngry PostFear
# Pre          1.00     -0.94        0.94    0.94     -0.74    -0.85
# PostHappy   -0.94      1.00       -0.93   -0.89      0.72     0.83
# PostNeutral  0.94     -0.93        1.00    0.87     -0.64    -0.78
# PostSad      0.94     -0.89        0.87    1.00     -0.80    -0.90
# PostAngry   -0.74      0.72       -0.64   -0.80      1.00     0.93
# PostFear    -0.85      0.83       -0.78   -0.90      0.93     1.00

# Note: Axis orientation in the final solution (on theta/delta) is arbitrary and correlations may be negative


# Export the delta correlations to LaTeX
R = abs(rcorr(locations)$r)
or = c("Pre","PostHappy","PostNeutral","PostFear","PostSad","PostAngry")
R = R[or,or]
cat(M.pre$toLatex(as.data.frame(R),caption="Correlations betwwen pre and post induction item peak locations across conditions.",label="tab:delta-corr", symmetric=TRUE))

# Emotion locations estimated accross conditions are very similar: ~88% of variance explained by a single component
pca1 = princomp(locations,cor=TRUE)

# Eigenvalues
e = pca1$sdev**2
e
#     Comp.1     Comp.2     Comp.3     Comp.4     Comp.5     Comp.6 
# 5.24437276 0.51205150 0.09858735 0.05991647 0.04613921 0.03893272 

# 87.4% of variance explained by a single component
e[1] / 6

#    Comp.1 
# 0.8794292 

#---------------------------------- Reanalysis with common emotion locations -------------------------------------

# Unfolding on separate Pre and Post data but with global locations (deltas) fixed
# This allows to get comparable subjects' locations between Pre- and Post-induction, and measure change

fixed.delta = M.pre$delta # Fixed deltas based on the merged data analysis

# Syntax: i) create the object, ii) provide fixed deltas as starting values, iii) declare them as fixed in the estimation call
M.post2.neutral = BUM$new(data=Post.sel[gNeutral,], inverted=sel.inv)
M.post2.neutral$delta = fixed.delta
M.post2.neutral$Estimate(fixed=c("delta"),display=FALSE)

M.post2.angry = BUM$new(data=Post.sel[gAngry,], inverted=sel.inv)
M.post2.angry$delta = fixed.delta
M.post2.angry$Estimate(fixed=c("delta"),display=FALSE)

M.post2.happy = BUM$new(data=Post.sel[gHappy,], inverted=sel.inv)
M.post2.happy$delta = fixed.delta
M.post2.happy$Estimate(fixed=c("delta"),display=FALSE)

M.post2.sad = BUM$new(data=Post.sel[gSad,], inverted=sel.inv)
M.post2.sad$delta = fixed.delta
M.post2.sad$Estimate(fixed=c("delta"),display=FALSE)

M.post2.fear = BUM$new(data=Post.sel[gFear,], inverted=sel.inv)
M.post2.fear$delta = fixed.delta
M.post2.fear$Estimate(fixed=c("delta"),display=FALSE)


#------------------------ Figure 8: Visualize the induction impact on subjects ----------------------------------
colors = c(blue="#4285f4", red="#db4437", yellow="#f4b400", green="#0f9d58", black="#333333")

library(mgcv)
library(dglm)

M.post.fixed = BUM$new(data=Post.sel, inverted=sel.inv)
M.post.fixed$delta = fixed.delta
M.post.fixed$Estimate(fixed="delta",display=FALSE)

subjects$Pre.score = M.pre$theta
subjects$Post.score = M.post.fixed$theta
write.csv2(subjects, file="subjects.csv", row.names=FALSE)

x11(height=10,width=7)

plot.condition = function(cond, gIdx, ...) {

  xc = M.pre$theta[gIdx]
  yc = M.post.fixed$theta[gIdx]
  cat(cond, "cor =", cor(xc, yc), "\n")
  plot(xc, yc, type="n", main=cond, xlab="Pre-induction valence", ylab="Post-induction valence", cex.main=1.3, cex.lab=1.2, ...)
  m = dglm(yc ~ xc, dformula = ~xc)
  tt = seq(min(xc), max(xc), len=100)
  mu = predict(m, newdata=data.frame(xc=tt))
  phi = exp(predict(m$dispersion.fit, newdata=data.frame(xc=tt)))
  polygon(c(tt, rev(tt)),
          c(mu - 1.96*sqrt(phi), rev(mu + 1.96*sqrt(phi))),
          col=gray(0.9), border=NA)
  abline(0, 1, lwd=2, col="gray", lty=2)
  segments(xc, xc, xc, yc, col="gray")
  lines(tt, mu, col=colors[1], lty=2, lwd=2)
  points(xc, yc, pch=19, cex=.7)
}

par(mfrow=c(3,2))
plot.condition("Joy",     gHappy,   xlim=c(-3,3), ylim=c(-4,4))
plot.condition("Neutral", gNeutral, xlim=c(-3,3), ylim=c(-4,4))
plot.condition("Sad",     gSad,     xlim=c(-3,3), ylim=c(-4,4))
plot.condition("Angry",   gAngry,   xlim=c(-3,3), ylim=c(-4,4))
plot.condition("Fear",    gFear,    xlim=c(-3,3), ylim=c(-4,4))

dev.copy2pdf(file="Images/Fig.8-induction-effect.pdf")


#------------------------ Figure 9: Reactivity of response curves to induction ---------------------------

th = seq(-2.5, 2.5, len=100)
items = c("Depressed","Anxious","Stressed","Bored","InControl","Active","Calm","Enthusiastic","Passive","Satisfied","Happy","Serene")
items=c("Depressed","Angry","Scared","Anxious","Stressed","Bored","InControl","Active","Calm","Enthusiastic","Joyful","Passive","Satisfied","Happy","Serene")
k = order(M.pre$delta[items])
items = items[k]

x11(height=10,width=7.5)
colors = c(blue="#4285f4", red="#db4437", yellow="#f4b400", green="#0f9d58", black="#333333")

fit.pre          = M.pre$Predict(tt=th)$mu
fit.post.happy   = M.post2.happy$Predict(tt=th)$mu
fit.post.neutral = M.post2.neutral$Predict(tt=th)$mu
fit.post.sad     = M.post2.sad$Predict(tt=th)$mu
fit.post.fear    = M.post2.fear$Predict(tt=th)$mu
fit.post.angry   = M.post2.angry$Predict(tt=th)$mu

par(mfrow=c(5,3), mar=c(4.5,4.5,2.5,2))
for(i in 1:length(items)) {

  y0 = fit.pre[,items[i]]
  y1 = fit.post.happy[,items[i]]
  y2 = fit.post.neutral[,items[i]]
  y3 = fit.post.sad[,items[i]]
  y4 = fit.post.fear[,items[i]]
  y5 = fit.post.angry[,items[i]]

  if(items[i] %in% sel.inv) {
    y0 = 1- y0
    y1 = 1- y1
    y2 = 1- y2
    y3 = 1- y3
    y4 = 1- y4
    y5 = 1- y5
  }

  plot(th,fit.pre[,items[i]],type="n",ylim=c(0,1.03),main=items[i],xlab=expression(theta),ylab="Response",cex.axis=1.1,cex.main=1.3,cex.lab=1.2)

  lines(th,y0,col="grey",lwd=5)
  lines(th,y1,col=colors[1],lwd=2,lty=2)
  lines(th,y2,col=colors[2],lwd=2,lty=3)
  lines(th,y3,col=colors[3],lwd=2,lty=4)
  lines(th,y4,col=colors[4],lwd=2,lty=5)
  lines(th,y5,col=colors[5],lwd=2,lty=6)

  if(i==1) legend("topright",bty="n",legend=c("Pre-ind.","Joy","Neutral","Sadness","Fear","Anger"),col=c("grey",colors),lwd=2,lty=1:6, cex=1)
}

dev.copy2pdf(file="Images/Fig.9-ICC-impact.pdf")

