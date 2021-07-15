# Orthogonal Matching Pursuit Method
The orthogonal matching pursuit (OMP) algorithm used in [this folder](https://github.com/jduffield65/iss/tree/master/%40iss_OMP_ConstantBackground_WeightDotProduct) differs from the [conventional algorithm](https://github.com/jduffield65/iss/blob/6b5cd2336e56ad844be8fe54cc36c38f8e0deba3/omp_free_background.m) in how the background is fitted and how the next atom (gene in this case) is selected. This document expalins how.

The function that carries out the OMP is [```o.get_omp_coefs```](https://github.com/jduffield65/iss/blob/6b5cd2336e56ad844be8fe54cc36c38f8e0deba3/@iss_OMP_ConstantBackground_WeightDotProduct/get_omp_coefs.m). It takes SpotColors [normalised by channel](https://github.com/jduffield65/iss/blob/6b5cd2336e56ad844be8fe54cc36c38f8e0deba3/@iss_Base/get_channel_norm.m) and fits the [background eigenvectors](https://github.com/jduffield65/iss/blob/6b5cd2336e56ad844be8fe54cc36c38f8e0deba3/@iss_OMP/get_background_codes.m) followed by successive genes, as long as they exceed a certain threshold.

## Fitting Background
The background is fitted using the function [```o.get_spot_residual_background```](https://github.com/jduffield65/iss/blob/6b5cd2336e56ad844be8fe54cc36c38f8e0deba3/@iss_OMP_ConstantBackground_WeightDotProduct/get_spot_residual_background.m). This assumes that there are ```o.nBP``` background eigenvectors and each one is just a strip in a color channel i.e. ```o.ompBackgroundChannelStrips=true```. An example for channel 1 is shown below:

<p float="left">
<img src="MathsImages/BackgroundCode.png" width = "450"> 
</p>

The background fitting procedure will be illustrated by fitting the above background vector to the following (channel normalised) spot color:

<p float="left">
<img src="MathsImages/SpotColor.png" width = "450"> 
</p>

The basic procedure, shown below, is to find a weight factor that normalises the contribution from each round so that no one round dominates. Multiply both background code and spot color by this weight. Then the final coefficient for the background code is the dot product of the weighted spot color with the weighted background code.

<p float="left">
<img src="MathsImages/BackgroundFit.png" width = "700"> 
</p>

The [weight factor](https://github.com/jduffield65/iss/blob/efa8542e0331e5337bdf05846755e1999b830b0d/%40iss_OMP_ConstantBackground_WeightDotProduct/get_spot_residual_background.m#L25-L30) for fitting the channel B background vector in round r is:

<img src="https://i.upmath.me/svg/W_%7BB%2Cr%7D%20%3D%20%5Cfrac%7B1%7D%7B(%7Cs_%7BB%2Cr%7D%7C%2B%5Clambda)%5E%5Csigma%7D" alt="W_{B,r} = \frac{1}{(|s_{B,r}|+\lambda)^\sigma}" />

<img src="https://i.upmath.me/svg/s_%7BB%2Cr%7D%7D" alt="s_{B,r}}" /> is the spot color in channel B, round r. <img src="https://i.upmath.me/svg/%5Clambda" alt="\lambda" /> is [```o.ompWeightShift```](https://github.com/jduffield65/iss/blob/0f7a804ea1b9f3d845788b826476ab459ba17388/%40iss_OMP_ConstantBackground_WeightDotProduct/iss_OMP_ConstantBackground_WeightDotProduct.m#L18) which is just a small number to stop <img src="https://i.upmath.me/svg/W_%7BB%2Cr%7D" alt="W_{B,r}" /> blowing up for small <img src="https://i.upmath.me/svg/s_%7BB%2Cr%7D%7D" alt="s_{B,r}}" />. <img src="https://i.upmath.me/svg/%5Csigma" alt="\sigma" /> is [```o.ompWeightPower```](https://github.com/jduffield65/iss/blob/0f7a804ea1b9f3d845788b826476ab459ba17388/%40iss_OMP_ConstantBackground_WeightDotProduct/iss_OMP_ConstantBackground_WeightDotProduct.m#L22), the lower <img src="https://i.upmath.me/svg/%5Csigma" alt="\sigma" />, the greater the contributions of the rounds where <img src="https://i.upmath.me/svg/s_%7BB%2Cr%7D%7D" alt="s_{B,r}}" /> is large. 

For <img src="https://i.upmath.me/svg/%5Clambda%3D0.01" alt="\lambda=0.01" /> and <img src="https://i.upmath.me/svg/%5Csigma%3D0.9" alt="\sigma=0.9" />, the weight factor is shown in the top row. The third row just shows the second row multiplied by the weight factor. If <img src="https://i.upmath.me/svg/%5Clambda%3D0" alt="\lambda=0" /> and <img src="https://i.upmath.me/svg/%5Csigma%3D1" alt="\sigma=1" />, then the weighted spot color shown would have an absolute value of 1 for all rounds in channel 1.

The [final coefficient](https://github.com/jduffield65/iss/blob/efa8542e0331e5337bdf05846755e1999b830b0d/%40iss_OMP_ConstantBackground_WeightDotProduct/get_spot_residual_background.m#L37) of the background vector for channel B is given by:
<img src="https://i.upmath.me/svg/C_B%20%3D%20%5Cfrac%7B%5Csum_%7Br%3D1%7D%5E7W_%7BB%2Cr%7D%5E2s_%7BB%2Cr%7Dg_%7BB%2Cr%7D%7D%7B%5Csum_%7Br%3D1%7D%5E7W_%7BB%2Cr%7D%5E2g_%7BB%2Cr%7D%5E2%7D" alt="C_B = \frac{\sum_{r=1}^7W_{B,r}^2s_{B,r}g_{B,r}}{\sum_{r=1}^7W_{B,r}^2g_{B,r}^2}" />

<img src="https://i.upmath.me/svg/g_%7BB%2Cr%7D" alt="g_{B,r}" /> is the value of the background vector for channel B in round r. In the procedure shown, this is equivalent to multiplying the two weighted codes in the third row together, then dividing by the square of the weighted background. The result of this is shown as the Weighted Dot Product in the bottom row. <img src="https://i.upmath.me/svg/C_B" alt="C_B" /> is then found by summing this. It is clear from the formula that if <img src="https://i.upmath.me/svg/s_%7BB%2Cr%7D%20%3D%20%5Cmu%20g_%7BB%2Cr%7D" alt="s_{B,r} = \mu g_{B,r}" /> then <img src="https://i.upmath.me/svg/C_B%20%3D%20%5Cmu" alt="C_B = \mu" /> for any value of <img src="https://i.upmath.me/svg/%5Cmu" alt="\mu" />.

The spot color found after fitting the background is shown in the middle plot below. The values of <img src="https://i.upmath.me/svg/C_b" alt="C_b" /> for all channels, b, are kept constant after this initial fitting which differs from the usual OMP method where they are re-fit after each subsequent gene is added. The reason for this is that sometimes an unusual background can be fit to justify further genes. 

<p float="left">
<img src="MathsImages/BackgroundResult.png" width = "700"> 
</p>

## Fitting Genes
After fitting the background, we need to decide which gene to fit next. This is done by the function [```o.get_weight_gene_dot_product2```](https://github.com/jduffield65/iss/blob/6b5e51b1bcc4f11ef7221cd2ffca18a2b45cbabf/@iss_OMP_ConstantBackground_WeightDotProduct/get_weight_gene_dot_product2.m). It finds a modified dot product <img src="https://i.upmath.me/svg/DP_g" alt="DP_g" /> for each gene, g, and then the next gene to be fit is the gene with the largest value. 

For the spot color shown above, the next gene to be fit is Aldoc, with a bled code shown below.

<p float="left">
<img src="MathsImages/AldocCode.png" width = "450"> 
</p>

The below plot illustrates the procedure for finding <img src="https://i.upmath.me/svg/DP_g" alt="DP_g" /> for Aldoc to be 6.098, with a particular focus on the contribution from round 7. The basic procedure is to find a dot product between the spot color and the gene code for each round and then summing these individual dot products with a weighting indicating the strength of the gene in each round.

<p float="left">
<img src="MathsImages/AldocFit.png" width = "1000"> 
</p>

First, for each gene, we need to find a weight for each round. This is based on the [gene efficiency](https://github.com/jduffield65/iss/blob/6b5e51b1bcc4f11ef7221cd2ffca18a2b45cbabf/@iss_OMP/get_gene_efficiencies.m) as shown in the top left plot above. The [formula](https://github.com/jduffield65/iss/blob/6b5e51b1bcc4f11ef7221cd2ffca18a2b45cbabf/%40iss_OMP_ConstantBackground_WeightDotProduct/get_omp_coefs.m#L34-L36) is:

<img src="https://i.upmath.me/svg/w_%7Bgr%7D%20%3D%20%5Csqrt%7B%5Cfrac%7B1%7D%7B1%2Be%5E%7B-%5Calpha(%5Cvarepsilon_%7Bg%2Cr%7D-%5Cbeta)%7D%7D%7D" alt="w_{gr} = \sqrt{\frac{1}{1+e^{-\alpha(\varepsilon_{g,r}-\beta)}}}" />

<img src="https://i.upmath.me/svg/%5Cvarepsilon_%7Bg%2Cr%7D" alt="\varepsilon_{g,r}" /> is the gene efficiency for gene g in round r. <img src="https://i.upmath.me/svg/%5Calpha" alt="\alpha" /> is [```ompNormBledCodeScale```](https://github.com/jduffield65/iss/blob/0f7a804ea1b9f3d845788b826476ab459ba17388/%40iss_OMP_ConstantBackground_WeightDotProduct/iss_OMP_ConstantBackground_WeightDotProduct.m#L32), the larger this, the sharper the step function. <img src="https://i.upmath.me/svg/%5Cbeta" alt="\beta" /> is [```ompNormBledCodeShift```](https://github.com/jduffield65/iss/blob/0f7a804ea1b9f3d845788b826476ab459ba17388/%40iss_OMP_ConstantBackground_WeightDotProduct/iss_OMP_ConstantBackground_WeightDotProduct.m#L31), this controls where the step function is centered. For the plot shown, <img src="https://i.upmath.me/svg/%5Calpha%3D7" alt="\alpha=7" /> and <img src="https://i.upmath.me/svg/%5Cbeta%3D0.5" alt="\beta=0.5" />.

We next square and [normalise this](https://github.com/jduffield65/iss/blob/6b5e51b1bcc4f11ef7221cd2ffca18a2b45cbabf/%40iss_OMP_ConstantBackground_WeightDotProduct/get_weight_gene_dot_product2.m#L14-L16) as shown in the top right plot. The formula is:

<img src="https://i.upmath.me/svg/W_%7Bgr%7D%5E2%20%3D%20%5Cfrac%7BnR%20%5Ctimes%20w_%7Bgr%7D%5E2%7D%7B%5Csum_%7Br%3D1%7D%5E%7BnR%7Dw_%7Bgr%7D%5E2%7D" alt="W_{gr}^2 = \frac{nR \times w_{gr}^2}{\sum_{r=1}^{nR}w_{gr}^2}" />

Here, nR = o.nRounds, and it is normalised such that <img src="https://i.upmath.me/svg/%5Csum_%7Br%3D1%7D%5E%7BnR%7DW_%7Bgr%7D%5E2%3DnR" alt="\sum_{r=1}^{nR}W_{gr}^2=nR" />.

Now we need to find the dot product for each round. For a particular round R, this is:

<img src="https://i.upmath.me/svg/DP_%7Bg%2CR%7D%20%3D%20%5Cfrac%7B%5Csum_%7Bb%3D1%7D%5E7s''_%7Bb%2CR%7Dg_%7Bb%2CR%7D%7D%7B%5CBigg(%5Csqrt%7B%5Csum_%7Bb%3D1%7D%5E7%7Bs''_%7Bb%2CR%7D%5E2%7D%7D%2B%5Clambda%5CBigg)%5E%7B%5Csigma%7D%5Csqrt%7B%5Csum_%7Bb%3D1%7D%5E7%7Bg_%7Bb%2CR%7D%5E2%7D%7D%7D" alt="DP_{g,R} = \frac{\sum_{b=1}^7s''_{b,R}g_{b,R}}{\Bigg(\sqrt{\sum_{b=1}^7{s''_{b,R}^2}}+\lambda\Bigg)^{\sigma}\sqrt{\sum_{b=1}^7{g_{b,R}^2}}}" />

Here, <img src="https://i.upmath.me/svg/s''" alt="s''" /> refers to the current spot residual, so for this example it is the spot color - background but if we are on the next iteration it would be spot color - background - first gene. We include both <img src="https://i.upmath.me/svg/%5Clambda" alt="\lambda" /> and <img src="https://i.upmath.me/svg/%5Csigma" alt="\sigma" /> again to stop a blow up and to give slightly greater contribution to the more intense rounds. If <img src="https://i.upmath.me/svg/%5Clambda%3D0" alt="\lambda=0" /> and <img src="https://i.upmath.me/svg/%5Csigma%3D1" alt="\sigma=1" />, then we are just finding the [cosine of the angle between the vectors](https://proofwiki.org/wiki/Cosine_Formula_for_Dot_Product).

In the procedure shown, the second row shows <img src="https://i.upmath.me/svg/s''_R" alt="s''_R" /> and <img src="https://i.upmath.me/svg/g_R" alt="g_R" /> with R=7. The third row shows the [spot weight](https://github.com/jduffield65/iss/blob/6b5e51b1bcc4f11ef7221cd2ffca18a2b45cbabf/%40iss_OMP_ConstantBackground_WeightDotProduct/get_weight_gene_dot_product2.m#L42-L54) which is:

<img src="https://i.upmath.me/svg/%5Cfrac%7B1%7D%7B%5CBigg(%5Csqrt%7B%5Csum_%7Bb%3D1%7D%5E7%7Bs''_%7Bb%2CR%7D%5E2%7D%7D%2B%5Clambda%5CBigg)%5E%7B%5Csigma%7D%7D" alt="\frac{1}{\Bigg(\sqrt{\sum_{b=1}^7{s''_{b,R}^2}}+\lambda\Bigg)^{\sigma}}" /> 

and the [gene weight](https://github.com/jduffield65/iss/blob/6b5e51b1bcc4f11ef7221cd2ffca18a2b45cbabf/%40iss_OMP_ConstantBackground_WeightDotProduct/get_weight_gene_dot_product2.m#L28-L40) which is:

<img src="https://i.upmath.me/svg/%5Cfrac%7B1%7D%7B%5Csqrt%7B%5Csum_%7Bb%3D1%7D%5E7%7Bg_%7Bb%2CR%7D%5E2%7D%7D%7D" alt="\frac{1}{\sqrt{\sum_{b=1}^7{g_{b,R}^2}}}" />.

The next row shows the two preceding rows multiplied by each other i.e.

<img src="https://i.upmath.me/svg/%5Cfrac%7Bs''_R%7D%7B%5CBigg(%5Csqrt%7B%5Csum_%7Bb%3D1%7D%5E7%7Bs''_%7Bb%2CR%7D%5E2%7D%7D%2B%5Clambda%5CBigg)%5E%7B%5Csigma%7D%7D" alt="\frac{s''_R}{\Bigg(\sqrt{\sum_{b=1}^7{s''_{b,R}^2}}+\lambda\Bigg)^{\sigma}}" /> and <img src="https://i.upmath.me/svg/%5Cfrac%7Bg_R%7D%7B%5Csqrt%7B%5Csum_%7Bb%3D1%7D%5E7%7Bg_%7Bb%2CR%7D%5E2%7D%7D%7D" alt="\frac{g_R}{\sqrt{\sum_{b=1}^7{g_{b,R}^2}}}" />.

The bottom left plot then shows these two multiplied by each other. Summing over all values in this image then gives <img src="https://i.upmath.me/svg/DP_%7Bg%2CR%7D%20%3D%200.894" alt="DP_{g,R} = 0.894" />. 

The final dot product is then given by:

<img src="https://i.upmath.me/svg/DP_g%20%3D%20%5Csum_%7Br%3D1%7D%5E7W_%7Bgr%7D%5E2DP_%7Bgr%7D" alt="DP_g = \sum_{r=1}^7W_{gr}^2DP_{gr}" />

In the bottom right plot shown, NormRoundWeight is <img src="https://i.upmath.me/svg/W_%7Bgr%7D%5E2" alt="W_{gr}^2" />, RoundDotProduct is <img src="https://i.upmath.me/svg/DP_%7Bgr%7D" alt="DP_{gr}" /> and RoundContribution is <img src="https://i.upmath.me/svg/W_%7Bgr%7D%5E2DP_%7Bgr%7D" alt="W_{gr}^2DP_{gr}" />.

## Finding coefficients for selected genes
Now we know the best gene to add, we need to find a coefficiet for it. This is done in the standard OMP least squares way. If we are considering the first gene to be added [then](https://github.com/jduffield65/iss/blob/d0892cb0c001b4f380e15443791aef9f5f26ad4e/%40iss_OMP_ConstantBackground_WeightDotProduct/get_spot_residual.m#L22):

<img src="https://i.upmath.me/svg/C_g%20%3D%20%5Cfrac%7B%5Csum_%7Bb%3D1%7D%5E7%5Csum_%7Br%3D1%7D%5E7s'_%7Bbr%7Dg_%7Bbr%7D%7D%7B%5Csum_%7Bb%3D1%7D%5E7%5Csum_%7Br%3D1%7D%5E7g_%7Bbr%7D%5E2%7D" alt="C_g = \frac{\sum_{b=1}^7\sum_{r=1}^7s'_{br}g_{br}}{\sum_{b=1}^7\sum_{r=1}^7g_{br}^2}" />

Where <img src="https://i.upmath.me/svg/s'" alt="s'" /> refers to the spot residual post background removal i.e. spot color - background. It is always this for every iteration as we refit the genes we have previously found when subsequent genes are added. 

If we are on a later iteration with more than one gene, then the MATLAB [backslash operation](https://github.com/jduffield65/iss/blob/d0892cb0c001b4f380e15443791aef9f5f26ad4e/%40iss_OMP_ConstantBackground_WeightDotProduct/get_spot_residual.m#L24) is used to solve:

<img src="https://i.upmath.me/svg/C_G%20%3D%20%5Cbigg(G%5ETG%5Cbigg)%5E%7B-1%7DG%5ETs'" alt="C_G = \bigg(G^TG\bigg)^{-1}G^Ts'" />

where G is a 49 x nGenesAdded matrix where each column is the bled code of a gene that is added and <img src="https://i.upmath.me/svg/C_G%20%3D%20%5BC_%7Bg1%7D%2C%20C_%7Bg2%7D%20...%5D" alt="C_G = [C_{g1}, C_{g2} ...]" />.

Genes are continually added in this way until there are [```o.ompMaxGenes```](https://github.com/jduffield65/iss/blob/d0892cb0c001b4f380e15443791aef9f5f26ad4e/%40iss_OMP/iss_OMP.m#L77-L79) added to explain a particular spot color or:

<img src="https://i.upmath.me/svg/%5CDelta%20Residual%20%3D%20Residual_%7Bi-1%7D%20-%20Residual_%7Bi%7D%20%3C%20ResidualThresh%20" alt="\Delta Residual = Residual_{i-1} - Residual_{i} &lt; ResidualThresh " />

Where:

<img src="https://i.upmath.me/svg/Residual%20%3D%20%5Csqrt%7B%5Csum_%7Bb%3D1%7D%5E7%5Csum_%7Br%3D1%7D%5E7s''_%7Bbr%7D%5E2%7D" alt="Residual = \sqrt{\sum_{b=1}^7\sum_{r=1}^7s''_{br}^2}" />

[ResidualThresh](https://github.com/jduffield65/iss/blob/6b5e51b1bcc4f11ef7221cd2ffca18a2b45cbabf/%40iss_OMP_ConstantBackground_WeightDotProduct/get_omp_coefs.m#L10-L17) is based on the second largest value in the initial spot color. The process is more easily explained through the example spot we have been using, this has ResidualThresh = 0.036:

<p float="left">
<img src="MathsImages/OMP_FullFit.png" width = "1000"> 
</p>

Here, i=1 refers to just Aldoc so <img src="https://i.upmath.me/svg/Residual_i%20%3D%200.46" alt="Residual_i = 0.46" />, <img src="https://i.upmath.me/svg/Residual_%7Bi-1%7D%20%3D%201.12" alt="Residual_{i-1} = 1.12" /> and <img src="https://i.upmath.me/svg/%5CDelta%20Residual%20%3D%200.66" alt="\Delta Residual = 0.66" />. <img src="https://i.upmath.me/svg/%5CDelta%20Residual%20%3E%20ResidualThresh" alt="\Delta Residual &gt; ResidualThresh" /> so we proceed and accept Aldoc. 

For i=2, <img src="https://i.upmath.me/svg/Residual_i%20%3D%200.24" alt="Residual_i = 0.24" />, <img src="https://i.upmath.me/svg/Residual_%7Bi-1%7D%20%3D%200.46" alt="Residual_{i-1} = 0.46" /> and <img src="https://i.upmath.me/svg/%5CDelta%20Residual%20%3D%200.22" alt="\Delta Residual = 0.22" />. <img src="https://i.upmath.me/svg/%5CDelta%20Residual%20%3E%20ResidualThresh" alt="\Delta Residual &gt; ResidualThresh" /> so we proceed and accept Trp53i11. 

For i=3, <img src="https://i.upmath.me/svg/Residual_i%20%3D%200.21" alt="Residual_i = 0.21" />, <img src="https://i.upmath.me/svg/Residual_%7Bi-1%7D%20%3D%200.24" alt="Residual_{i-1} = 0.24" /> and <img src="https://i.upmath.me/svg/%5CDelta%20Residual%20%3D%200.03" alt="\Delta Residual = 0.03" />. <img src="https://i.upmath.me/svg/%5CDelta%20Residual%20%3C%20ResidualThresh" alt="\Delta Residual &lt; ResidualThresh" /> so we reject Hapln1 and end.

