TGraphAsymmErrors *slidingWindowAverage(TGraphAsymmErrors *input, int slidingWindow) {
    TGraphAsymmErrors *out  = (TGraphAsymmErrors *) input->Clone();
    bool isLogPlot = TString(input->GetName()).Contains("smcls");
    for (int i = 0, n = input->GetN(); i < n; ++i) {
        /*if (i >= 1 && i < n-1 && (input->GetX()[i] < input->GetX()[i-1] || input->GetX()[i] > input->GetX()[i+1])) {
            out->GetX()[i] = 0.5*(input->GetX()[i-1] + input->GetX()[i+1]);
            out->GetY()[i] = 0.5*(input->GetY()[i-1] + input->GetY()[i+1]);
            out->GetEYlow()[i] = 0.5*(input->GetEYlow()[i-1] + input->GetEYlow()[i+1]);
            out->GetEYhigh()[i] = 0.5*(input->GetEYhigh()[i-1] + input->GetEYhigh()[i+1]);
            continue;
        }*/
        double y0 = input->GetY()[i];
        double sum = 0, sumhi = 0, sumlo = 0, sumw = 0;
        for (int j = i-slidingWindow; j <= i+slidingWindow; ++j) {
            if (j < 0 || j >= n) continue;
            double y = input->GetY()[j], w = 1.0; // /(1.0 + abs(i-j));
            if (isLogPlot) {
                if (fabs(log(y0/y)) > log(2)) continue;    
            } else {
                if (fabs(y0-y) > 0.1*y0) continue;    
            }
            if (isLogPlot) {
                sum   += w*log(y);
                sumlo += w*log(input->GetEYlow()[j]);
                sumhi += w*log(input->GetEYhigh()[j]);
            } else {
                sum   += w*y;
                sumlo += w*input->GetEYlow()[j];
                sumhi += w*input->GetEYhigh()[j];
            }
            sumw  += w;
        }
        if (sumw == 0) continue;
        if (isLogPlot) {
            out->GetY()[i] = exp(sum/sumw);
            out->GetEYlow()[i] = exp(sumlo/sumw);
            out->GetEYhigh()[i] = exp(sumhi/sumw);
        } else {
            out->GetY()[i] = sum/sumw;
            out->GetEYlow()[i] = sumlo/sumw;
            out->GetEYhigh()[i] = sumhi/sumw;
        }
    }
    return out;
}

TGraphAsymmErrors *smoothSMCLs(TGraphAsymmErrors *input, int slidingWindow, int order) {
    double xi[99], yi[99];
    TString who = TString(input->GetName());
    bool nonzero   = true; //TString(input->GetName()).Contains("cls");
    bool isLogPlot = (who.Contains("smcls") || who.Contains("pval"));
    bool absolute  = (who.Contains("smcls") || who.Contains("pval"));
    double maxdel  = (who.Contains("smcls") || who.Contains("pval")) ? 0 : log(1.25); 
    TGraphAsymmErrors *inp = (TGraphAsymmErrors *) input->Clone();
    TGraphAsymmErrors *out = (TGraphAsymmErrors *) input->Clone();
    if (absolute) {
        for (int i = 0, n = input->GetN(); i < n; ++i) {
            inp->GetEYlow() [i] = inp->GetY()[i] - inp->GetEYlow()[i];
            inp->GetEYhigh()[i] = inp->GetY()[i] + inp->GetEYhigh()[i];
            out->GetEYlow() [i] = out->GetY()[i] - out->GetEYlow()[i];
            out->GetEYhigh()[i] = out->GetY()[i] + out->GetEYhigh()[i];
        }
    } 
    for (int which = -1; which <= +1; ++which) {
        double *inY = (which == 0 ? inp->GetY() : (which > 0 ? inp->GetEYhigh() : inp->GetEYlow()));
        double *ouY = (which == 0 ? out->GetY() : (which > 0 ? out->GetEYhigh() : out->GetEYlow()));
        for (int i = 0, n = input->GetN(); i < n; ++i) {
            double y0 = inY[i], x0 = input->GetX()[i]; int points = 0;
            if (nonzero && y0 == 0) continue;
            for (int j = i-slidingWindow; j <= i+slidingWindow; ++j) {
                if (j < 0 || j >= n) continue;
                if (nonzero && inY[j] == 0) continue;
                if (maxdel > 0) {
                    double del = inp->GetY()[j]/inp->GetY()[i]; 
                    if (fabs(log(del)) > maxdel) continue;
                }
                xi[points] = input->GetX()[j];
                yi[points] = isLogPlot ? log(inY[j]/y0) : inY[j];
                points++;
            }
            double ynew = isLogPlot ? 1 : y0;
            if (points > order+2) {
                ynew = smoothWithPolyFit(x0, order+1, points, xi, yi);
            } else if (points > 2) {
                ynew = smoothWithPolyFit(x0, 1, points, xi, yi);
            } else if (maxdel == 0 && nonzero && points == 1) {
                ouY[i] = 0; continue; // kill the blip
            } else continue; // nothing to do
            ouY[i] = isLogPlot ? y0 * exp(ynew) : ynew;
        }
    }
    if (absolute) {
        for (int i = 0, n = input->GetN(); i < n; ++i) {
            out->GetEYlow() [i] = +out->GetY()[i] - out->GetEYlow()[i];
            out->GetEYhigh()[i] = -out->GetY()[i] + out->GetEYhigh()[i];
        }
    }
    delete inp;
    return out;
}

TGraphAsymmErrors * removeGlitches(TGraphAsymmErrors *out) {
    do {
        int bad = -1; bool hasgood = false;
        for (int i = 0; i < out->GetN(); ++i) {
            if (out->GetEYlow()[i] == 0 && out->GetEYhigh()[i] == 0) {
                bad = i;
            } else { 
                hasgood = true; 
            }
        }
        if (!hasgood) return out;
        if (bad == -1) return out;
        out->RemovePoint(bad); 
    } while (1);
}

int nmasses = 0;
double masses[200];
void loadMasses(const char *file="masses.txt") {
    FILE *f = fopen(file,"r");
    float mass;
    while (fscanf(f,"%f", &mass) == 1) {
        masses[nmasses++] = mass;
    }
}

TGraphAsymmErrors *missingPoints(TGraphAsymmErrors *points) {
    if (!track_missing) return 0;
    if (points == 0) return 0;
    int n = points->GetN(); Double_t *xj = points->GetX(); Double_t *yj = points->GetY();
    double xmin = xj[0], xmax = xj[n-1];
    bool logint = true;
    TGraphAsymmErrors *ret = new TGraphAsymmErrors(); ret->SetName(TString("missing_")+points->GetName());
    // check if we have half-integer points
    bool halfint = false;
    for (int i = 0; i < n; ++i) { 
        if (points->GetX()[i] - floor(points->GetX()[i]) > 0.4) { halfint = true; break; }
    } 
    int nmiss = 0;
    for (int i = 0; i < nmasses; ++i) {
        double x = masses[i];
        if (masses[i] < xmin || masses[i] > xmax) continue;
        if ((x - floor(x)) > 0.4 && !halfint) continue;
        bool found = false; double xlo = 0, ylo = 0, xhi = 0, yhi = 0;
        for (int j = 0; j < n; ++j) {
            if (xj[j] < x) { xlo = xj[j]; ylo = yj[j]; }
            else if (xj[j] == x) { found = true; break; }
            else if (xhi == 0) { xhi = xj[j]; yhi = yj[j]; break; }
        }
        if (!found) {
            double y = ( yhi * (x-xlo) + ylo * (xhi - x) ) / (xhi - xlo);
            if (yhi > 0 && ylo > 0 && logint) {
                y = exp( ( log(yhi) * (x-xlo) + log(ylo) * (xhi - x) ) / (xhi - xlo) );
            }
            //printf("Interpolated missing point %3d from (%4.0f, %7.3f) + (%4.0f, %7.3f) --> (%4.0f, %7.3f)\n", masses[i], xlo, ylo, xhi, yhi, x, y);
            nmiss++;
            ret->Set(nmiss);
            ret->SetPoint(nmiss-1, x, y);
        }
    }
    ret->SetMarkerStyle(20);
    ret->SetMarkerSize(0.7);
    ret->SetMarkerColor(100);
    if (nmiss == 0) { delete ret; return 0; }
    return ret;
}

void findCrossings(TString who, TString xname, double threshold, double xmin, double xmax, TString what="95% CL limit on #sigma/#sigma_{SM}") {
    TGraphAsymmErrors *obs = (TGraphAsymmErrors *) gFile->Get(who);
    if (who == 0) return;
    int ihigh = -1, ilow = -1;
    for (int i = 0, n = obs->GetN(); i < n; ++i) {
        if (xmin <= obs->GetX()[i] && obs->GetX()[i] <= xmax) {
            if (xname.Contains("low")) {
                if (obs->GetY()[i] < threshold) {
                    ilow = i;
                    break;
                } else { 
                    ihigh = i;
                }
            } else  {
                if (obs->GetY()[i] < threshold) {
                    ilow = i;
                } else { 
                    ihigh = i;
                    break;
                }
            }
        }
    }
    if (ilow == -1 || ihigh == -1) { std::cout << "didn't find points." << std::endl; return; }
    double x1 = obs->GetX()[ilow], x2 = obs->GetX()[ihigh];
    double y1 = obs->GetY()[ilow], y2 = obs->GetY()[ihigh];
    TF1 *linear  = new TF1("linear", "[0] * (x - [1]) + [2]", xmin, xmax); 
    linear->SetParameters((y2-y1)/(x2-x1), x1, y1);
    TF1 *linlog = new TF1("linlog", "[0] * pow([1], (x - [2])/[3])", xmin, xmax);
    linlog->SetParameters(y1, y2/y1, x1, x2-x1);
    TF1 *loglog  = new TF1("linear", "[0] * pow([1], log(x/[2])/[3])", xmin, xmax);
    loglog->SetParameters(y1, y2/y1, x1, log(x2/x1));
    loglog->SetLineWidth(3); loglog->SetLineColor( 63); 
    linlog->SetLineWidth(5); linlog->SetLineColor(210); linlog->SetLineStyle(2);
    linear->SetLineWidth(3); linear->SetLineColor(  1);  linear->SetLineStyle(9);
    obs->SetMarkerStyle(21); obs->SetMarkerSize(1.3);
    double x_linear = linear->GetX(threshold, xmin, xmax);
    double x_linlog = linlog->GetX(threshold, xmin, xmax);
    double x_loglog = loglog->GetX(threshold, xmin, xmax);
    TH1D frame0("frame","frame", 1, xmin, xmax); frame0.Draw(); gStyle->SetOptStat(0);
    linear->Draw("SAME");  linlog->Draw("SAME");     loglog->Draw("SAME");
    obs->Draw("PX");    
    double ymin = threshold/2, ymax = threshold * 3;
    double xleg = 0.67, yleg = 0.18;
    if (who.Contains("smcls")) { 
        ymax = 3*threshold; ymin = 0;  
        xleg = 0.27; yleg = 0.55;
    }
    if (xname.Contains("low")) {
        xleg = 0.67 ;
        yleg = 0.66;
    }
    leg = newLegend(xleg,yleg,xleg+.28,yleg+.25); leg->SetTextSize(0.04);
    leg->AddEntry(linear, Form("Linear,  m_{H} = %.1f", x_linear), "L");
    leg->AddEntry(linlog, Form("Lin-log, m_{H} = %.1f", x_linlog), "L");
    leg->AddEntry(loglog, Form("Log-log, m_{H} = %.1f", x_loglog), "L");
    setCanvas(&frame0, "", ymin, ymax, what);
    TString myspam = "#splitline{"+SPAM+"}{"+channelFromName(who)+"}";
    finalize(who+"_crossing_"+xname,xmin,xmax,ymin,ymax, myspam);
}
 

TGraphAsymmErrors *getNSigmaLine(TGraphAsymmErrors *g, int sigma/* = +/-1 */) {
    if (g == 0) return 0;
    int n = g->GetN();
    TGraphAsymmErrors *ret = new TGraphAsymmErrors(n);
    for (int i = 0; i < n; ++i) {
        ret->SetPoint(i, g->GetX()[i], g->GetY()[i] + (sigma > 0 ? g->GetErrorYhigh(i) : - g->GetErrorYlow(i)));
        ret->SetPointError(i, 0, 0, 0, 0);
    }
    ret->SetLineWidth(g->GetLineWidth() > 1 ? g->GetLineWidth() - 1 : 1);
    ret->SetLineColor(g->GetLineColor());
    ret->SetLineStyle(2);
    return ret;
}


