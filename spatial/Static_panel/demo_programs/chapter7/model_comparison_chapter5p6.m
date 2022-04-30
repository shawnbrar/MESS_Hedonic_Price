% model_comparison_chapter5p6
clear all;

fnames = ['VERIZON COMMUNICATIONS INC  '
    'CENTURYTEL INC              '
    'CINCINNATI BELL INC         '
    'COMCAST CORP                '
    'DISNEY (WALT) CO            '
    'GRAY TELEVISION INC         '
    'HICKORY TECH CORP           '
    'TELEPHONE & DATA SYSTEMS INC'
    'SPRINT NEXTEL CORP          '
    'CABLEVISION SYS CORP  -CL A '
    'GENERAL COMMUNICATION  -CL A'
    'CBS CORP                    '
    'US CELLULAR CORP            '
    'EW SCRIPPS  -CL A           '
    'EMMIS COMMUNICTNS CP  -CL A '
    'NTN BUZZTIME INC            '
    'ATLANTIC TELE-NETWORK INC   '
    'USA MOBILITY INC            '
    'TIME WARNER INC             '
    'SAGA COMMUNICATIONS  -CL A  '
    'SINCLAIR BROADCAST GP  -CL A'
    'ECHOSTAR COMMUN CORP  -CL A '
    'PREMIERE GLOBAL SERVICES INC'
    'IDT CORP                    '
    'KVH INDUSTRIES INC          '];

firm_names = strvcat(fnames);

data = load('../demo_data/industry32.data');
% variables are
% 1 gvkey
% 2 year
% 3 ff_48
% 4 mkt leverage
% 5 log(sales)
% 6 roa
% 7 mtb
% 8 tangibility
% 9 zipcode area
% 10 latt of firm headquarters
% 11 long  of firm headquarters
% 12 countyfips
% 13 statefips

    T = 18; %years
    num_firms = size(data,1)/T;
    N = num_firms;
    nvar = 4;

% row-normalized Hoberg-Philips W-matrix of product market competitors
% for SIC 32 firms    
W_hp = load('../demo_data/Windustry32.data');
% create a binary matrix of competitor firms
Wbinary = zeros(N,N);
for i=1:N
    Wbinary(:,i) = (W_hp(i,:) > 0);
end
W_hp_equal = normw(Wbinary);
    
    % rearrange the data
    % all firms for year 1
    % all firms for year 2
    % etc
%     id = data(:,1);
%     idu = unique(id);
%     nu = length(idu);
%     index = data(:,1);
    year = data(:,2);
    yru = unique(year);
    nyr = length(yru);
    
    time_names = strvcat(num2str(yru));

    new_data = [];
    for i=1:nyr
        yri = yru(i,1);
        ind = find(year == yri);
        new_data = [new_data
                   data(ind,:)];
               
    end
    
    id = new_data(:,1);
    idu = unique(id);
    
% 3 Fama-French_48
% 4 mkt leverage
% 5 log(sales)
% 6 roa
% 7 mtb
% 8 tangibility
% 9 zipcode area
% 10 latt
% 11 long

mkt_leverage = new_data(:,4);
sales = new_data(:,5);
roa = new_data(:,6);
mtb = new_data(:,7);
tangibility = new_data(:,8);

    y = mkt_leverage; % return on assets
    x = [roa sales  mtb tangibility];
    
% ===================================================
model = 1; % model with firm and effects
[ywith,xwith,meanny,meannx,meanty,meantx]=demean(y,x,N,T,model);

info.lflag = 0;
result1 = lmarginal_static_panel(ywith,xwith,W_hp,N,T,info); 
in.width = 10000;
in.fmt = '%10.4f';
out = [result1.lmarginal];

result2 = lmarginal_static_panel(ywith,xwith,W_hp_equal,N,T,info); 
out = [out
       result2.lmarginal];
   
mprobs = model_probs(out);


% reshape for pretty print
prt_out = reshape(out,3,2);
prt_probs = reshape(mprobs,3,2);

out = [prt_out(:,1) prt_probs(:,1) prt_out(:,2) prt_probs(:,2)];


fprintf(1,'probs for models based  firm effects \n');
in.cnames = strvcat('Whp lmarginal','Whp prob' ,'Whp equal lmarginal','Whp equal prob');
in.rnames = strvcat('Model','slx','sdm','sdem');
in.width = 10000;
in.fmt = '%10.4f';
mprint(out,in);


prior.model = 1;
prior.rval = 5;
ndraw = 2500;
nomit = 500;

result1 = sdem_panel_FE_g(y,x,W_hp,T,ndraw,nomit,prior);
vnames = strvcat('y=mkt_leverage','x1=roa','x2=sales','x3=market2book','x4=tangibility');
prt_panel(result1,vnames);

result2 = sdem_panel_FE_g(y,x,W_hp_equal,T,ndraw,nomit,prior);
vnames = strvcat('y=mkt_leverage','x1=roa','x2=sales','x3=market2book','x4=tangibility');
prt_panel(result2,vnames);

result3 = sdm_panel_FE_g(y,x,W_hp,T,ndraw,nomit,prior);
vnames = strvcat('y=mkt_leverage','x1=roa','x2=sales','x3=market2book','x4=tangibility');
prt_panel(result3,vnames);


