

%%%% INPUT %%%
m = 10;                                         % # vnitřních bodů
overlap = 1;                                    % velikost překryvu


h = 2/(m);                                      % velikost kroku
f = -5 * ones(m-1, 1);
x = linspace(-1, 1, m + 1);
u_history = zeros(m+1,m);                       % uložení řešení z \Omega_1 na \Omega
v_history = zeros(m+1,m);                       % uložení řešení z \Omega_2 na \Omega
e1 = zeros(m+1,m);                              % uložení chyby na \Omega_1
e2 = zeros(m+1,m);                              % uložení chyby na \Omega_2


%definice tridiagonální matice
A = (1 / h^2) * diag(-2 * ones(m-1, 1)) + ...   %vektor o velikosti m-1
    (1 / h^2) * diag(ones(m-2, 1), 1) + ...     %vektor o velikosti m-2 (bez diagonály) - naddiagonála
    (1 / h^2) * diag(ones(m-2, 1), -1);         %vektor o velikosti m-2 (bez diagonály) - poddiagonála

t = m-1-overlap;                                % # bodů mimo overlap
m_Omg1 = overlap+t/2;                            % konec \Omega_1
m_pg = overlap;
m_mg = m_Omg1 - overlap +1;


A1 = A(1:m_Omg1,1:m_Omg1);   
A2 = A(m_mg:end,m_mg:end);

u = zeros(m_Omg1,1);
v = zeros(m_Omg1,1);
w = zeros(1,m_Omg1);
w(end) = 1;
r = zeros(1,m_Omg1);
r(1)=1;

for k = 1:m
    u = A1\([f(1:m_Omg1-1); f(m_Omg1) - 1/(h^2)*(v(m_pg+1))]);
    u_full = [0; u; v(m_pg+1); nan(t/2, 1)];
    v = A2\( [f(m_mg) - 1/(h^2)*(u(m_mg-1));f(m_mg+1:m-1)]);
    v_full = [nan(t/2, 1); u(m_mg-1); v; 0];

    u_history(:, k) = u_full;
    v_history(:, k) = v_full;
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
