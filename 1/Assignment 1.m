

% 3- Test the quantizer/dequantizer functions
n_bits = 3;
xmax = 6;
x = -6:0.01:6;
m_values = [0, 1];
for m = m_values
    q_ind = UniformQuantizer(x, n_bits, xmax, m);
    deq_val = UniformDequantizer(q_ind, n_bits, xmax, m);
    figure;
    plot(x, x, 'b', x, deq_val, 'r', x, q_ind, 'g');
    title(sprintf('m=%d', m));
    legend('Original Signal', 'Quantized-Dequantized Signal', 'Quantization Indices');
end

% 4- Test on random input signal
sim_snr=[];
thr_snr=[];
x = random('Uniform',-5,5,1,10000);
xmax=max(abs(x));
m=0;
for n_bits= 2:1:8
  y=UniformQuantizer(x,n_bits,xmax,m);
  y_deq = UniformDequantizer(y, n_bits, xmax, m);
  error=abs(x-y_deq);
  sim_snr = [sim_snr, mean(x.^2)/mean(error.^2)];
  scale=(3*((2^n_bits)^2))/(xmax^2);
  thr_snr = [thr_snr, scale*mean(x.^2)];
end
n_bits= 2:1:8;


figure;
plot(n_bits, mag2db(sim_snr),'b', n_bits,mag2db(thr_snr),'r');
xlabel('Number of bits');
ylabel('SNR (dB)');
title('SNR_theo and  SNR_sim Number of Bits');
legend('sim','thr');

% 5- Test on non-uniform random input
sim_snr=[];
thr_snr=[];
size = [1 10000];
x_exp = exprnd(1,size);
sign = (randi([0,1],size)*2)-1;
x = x_exp.*sign;
xmax=max(abs(x));
m=0;
for n_bits= 2:1:8
y=UniformQuantizer(x,n_bits,xmax,m);
y_deq = UniformDequantizer(y, n_bits, xmax, m);
error=abs(x-y_deq);
sim_snr = [sim_snr, mean(x.^2)/mean(error.^2)];
scale=(3*((2^n_bits)^2))/(xmax^2);
thr_snr = [thr_snr, scale*mean(x.^2)];
end
n_bits= 2:1:8;
figure;
plot(n_bits, mag2db(sim_snr),'b', n_bits,mag2db(thr_snr),'r');
xlabel('Number of bits');
ylabel('SNR (dB)');
title('SNR_theo and  SNR_sim Number of Bits');
legend('sim','thr');
% 6- Non-uniform Mu-law quantization

figure();
x_norm=x/xmax;
c={"r",'g','b','m'};
i=1;
for mu=[0, 5, 100,200]
sim_snr=[];
thr_snr=[];
if(mu~=0)
x_comp = Compression(x_norm,mu,sign);
else
x_comp=x;
end
ymax=max(abs(x_comp));
for n_bits= 2:1:8
y = UniformQuantizer(x_comp, n_bits,ymax, m);
y_deq = UniformDequantizer(y, n_bits, ymax, m);
if(mu~=0)
y_expand = Expansion(y_deq,mu,sign);
y_deq = y_expand *xmax;
end
error=abs(x-y_deq);
sim_snr = [sim_snr, mean(x.^2)/mean(error.^2)];
if(mu~=0)
scale=(3*((2^n_bits)^2));
thr_snr = [thr_snr,scale/((log(1+mu))^2)];
else
scale=(3*((2^n_bits)^2))/(xmax^2);
thr_snr = [thr_snr, scale*mean(x.^2)];
end
end
n_bits= 2:1:8;
% plotTwoOutput(n_bits,mag2db(sim_snr),n_bits,mag2db(thr_snr),strcat('SNR(6) mu= ',
% num2str(mu)),"n-bits","snr");
% legend('sim','thr');
plot(n_bits,mag2db(sim_snr),'-','color',c{i})
hold on
plot(n_bits,mag2db(thr_snr),'--','color',c{i})
title("SNR (6)")
xlabel("n-bits")
ylabel("snr")
i=i+1;
end
legend('sim mu=0','thr mu=0','sim mu=5','thr mu=5','sim mu=100','thr mu=100','sim mu=200','thrmu=200');
function y = Compression(x, u, sign)
y=sign .* (log(1+u*abs(x))/log(1+u));
end
function y = Expansion(x, u, sign)
y= sign .*(((1+u).^abs(x)-1)/u);
end
function q_ind = UniformQuantizer(in_val, n_bits, xmax, m)
levels = 2 ^ n_bits;
delta = 2 * xmax / levels;
q_ind = floor((in_val - ((m) * (delta / 2) - xmax)) / delta);
q_ind(q_ind<0) = 0;
end

function deq_val = UniformDequantizer(q_ind, n_bits, xmax, m)
levels = 2 ^ n_bits;
delta = 2 * xmax / levels;
deq_val = ((q_ind) * delta) + ((m+1) * (delta / 2) - xmax);
end
