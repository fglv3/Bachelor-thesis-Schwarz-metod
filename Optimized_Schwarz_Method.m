
%%% INPUT %%%
m = 10;                                        % počet vnitřních bodů + pravý okraj
overlap = 1;                                   % velikost překryvu

f = -3 * ones(m-1, 1);
h = 2 / m;                                     % velikost kroku
x = linspace(-1, 1, m + 1);
u_history = zeros(m+1,m);                      % uložení řešení z \Omega_1 na \Omega
v_history = zeros(m+1,m);                      % uložení řešení z \Omega_2 na \Omega
e1 = zeros(m+1,m);                             % uložení chyby na \Omega_1
e2 = zeros(m+1,m);                             % uložení chyby na \Omega_1

% Matice A
A = (1 / h^2) * diag(-2 * ones(m-1, 1)) + ...
    (1 / h^2) * diag(ones(m-2, 1), 1) + ...
    (1 / h^2) * diag(ones(m-2, 1), -1);

t = m - 1 - overlap;                          % počet bodů mimo overlap
m_Omg1 = overlap + t/2 + 1;                   % # neznávým v Omg1
m_pg = overlap + 2;                           % index konce Omg1
m_mg = m_Omg1-overlap-1;                      % index začátek Omg2

% Přesné řešení
u_ex = [0; A\f; 0];

w = zeros(1,m_Omg1); 
w(end) = 1;
r = zeros(1,m_Omg1); 
r(1) = 1;

% Volba parametru 
p = 1 / (1-(overlap+1) * h/2);

% rozdělení matice A dle podmnožin Omwga_i s Robinovými podmínkami
A1 = A(1:m_Omg1, 1:m_Omg1);
A1(end,:) = 0;
A1(end,end) =1/h+p;
A1(end,end-1) = -1/h;

A2 = A(m_mg:end, m_mg:end);
A2(1,:) = 0;
A2(1,1) = (p+1/h);
A2(1,2) = -1/h;

v = zeros(m_Omg1,1);                        % počáteční aproximace

for k = 1:m
    robin_1 = (v(m_pg+1)-v(m_pg-1) )/ (2*h) + p * v(m_pg);
    u = A1 \ ([f(1:m_Omg1-1); robin_1-f(m_Omg1)*(h/2)]);

    robin_2 = (u(m_mg-1)- u(m_mg+1)) / (2*h) + p * u(m_mg);
    v = A2 \ ([robin_2 - f(m_mg)*(h/2); f(m_mg+1:m-1)]);

    u_full = [0; u; nan(t/2,1)];
    v_full = [nan(t/2,1); v; 0];

    u_history(:,k) = u_full;
    v_history(:,k) = v_full;

end


figure
hold on
% vykreslení všech iterací Schwarzovy metody pro obě poddomény
for k = 1:m
    plot(x, u_history(:,k), 'b', 'HandleVisibility', 'off');
    plot(x, v_history(:,k), 'r', 'HandleVisibility', 'off');
end

% Přidání zástupců do legendy pro Schwarz iterace
plot(nan, nan, 'b', 'DisplayName', '$u_1^n$');
plot(nan, nan, 'r', 'DisplayName', '$u_2^n$');
legend('Interpreter', 'latex', 'Location', 'best');
grid on;

legend("Position", [0.77147,0.79372,0.11429,0.10595])

