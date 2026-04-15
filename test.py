import numpy as np
import sympy as sp
import mesh_2d as m2d
import system_2d as s2d
import finite_element_2d as f2e
import graph2d as g2d



# # test1
# def ft1(x, y, A=100, x0=0.5, y0=0.5, const=0.1):
#     return A * np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / const)


# def ft2(x, y, A=100, x0=0.5, y0=0.5, const=1):
#     return A * np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / const)


# test2
# sourcest3 = [(0.5, 0.5)]
# strengthst3 = [100]

# sourcest4 = [(0.5, 0.5), (0.2, 0.8)]
# strengthst4 = [100, 50]

# atol = 0.05


# def ft3(x, y): return sum(
#     s for (x0, y0), s in zip(sourcest3, strengthst3)
#     if np.isclose(x, x0, atol=atol) and np.isclose(y, y0, atol=atol)
# )


# def ft4(x, y): return sum(
#     s for (x0, y0), s in zip(sourcest4, strengthst4)
#     if np.isclose(x, x0, atol=atol) and np.isclose(y, y0, atol=atol)
# )


# test3
# def k1t5(x, y):
#     return 0.1 if 0.45 < x < 1.55 else 10.0


# def k2t5(x, y):
#     return 1.0


# def ft5(x, y):
#     return 100 * np.exp(-((x - 0.5) ** 2 + (y - 0.5) ** 2) / 0.001)


# def k1t6(x, y):
#     if 0.45 < x < 1.55:
#         return 0.05
#     return 1.0


# def k2t6(x, y):
#     if 0.65 < y < 1.75:
#         return 0.01
#     return 1.0


# def ft6(x, y):
#     if x < 0.3 and y < 0.3:
#         return 20
#     elif x > 0.7 and y > 0.7:
#         return 5
#     else:
#         return 0


# ver
# def fv(x, y):
#     return 2 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)
# def fv(x, y):
#     return 2 * sp.pi ** 2 * sp.sin(sp.pi * x) * sp.sin(sp.pi * y)

# Нове f для реакції-дифузії з beta = const
# du/dt - Δu - beta*u = f
# => f = (-1 + 2π² - beta) * e^{-t} * sin(πx) * sin(πy)

# def exact_solution(x, y):
#     return np.sin(np.pi * x) * np.sin(np.pi * y)


def ug_1(x, y):
    return 0
    # return y**2
    # return y
    # return y**2


def ug_2(x, y):
    return 0
    # return y**2
    # return y+1
    # return y**2 + 1


def ug_3(x, y):
    return 0
    # return x**2
    # return x


def ug_4(x, y):
    return 0
    # return x**2+1
    # return x+1


def k1(x, y):
    return 1


def k2(x, y):
    return 1

def beta(x, y):
    return 0.5

def fv_space_reaction(x, y, beta_val=0.5):
    return (2 * np.pi**2 - 1 + beta_val) * np.sin(np.pi * x) * np.sin(np.pi * y)

def fv_space(x, y):
    """Просторова (базова) частина функції джерела"""
    return (2 * np.pi**2 - 1) * np.sin(np.pi * x) * np.sin(np.pi * y)

def exact_solution_t(x, y, t):
    """Точний аналітичний розв'язок, що залежить від часу"""
    return np.exp(-t) * np.sin(np.pi * x) * np.sin(np.pi * y)


def main():
    # verticest1 = [(0, 0), (1.5, 0), (1, 1), (0, 0.75)]
    # verticest2 = [(0, 0), (1.5, 0), (1, 1), (0, 0.75)]
    # verticest3 = [(0, 0), (1, 0), (1.1, 1), (0, 0.9)]
    # verticest4 = [(0, 0), (1, 0), (1.1, 1), (0, 0.9)]
    # verticest5 = [(0, 0), (1.5, 0), (1, 1), (0, 0.75)]
    # verticest6 = [(0, 0), (1, 0), (1, 1), (0, 1)]
    verticesv = [(0, 0), (1, 0), (1, 1), (0, 1)]

    ap = 3

    ug = [
        [f2e.TypeOfBoundCond.DIRICHLET, ug_3],  # нижнє (y = 0)
        [f2e.TypeOfBoundCond.DIRICHLET, ug_2],  # ліве (x = 0)
        [f2e.TypeOfBoundCond.DIRICHLET, ug_4],  # праве (x = 1)
        [f2e.TypeOfBoundCond.DIRICHLET, ug_1]  # верхнє (y = 1)
    ]

    element_type = 'D2QU4N'
    p = 10
    m = 10

    # ver
    print("\n--- Запуск верифікаційного тесту ---")
    NLv, ELv = m2d.uniform_mesh_with_vertices(verticesv, p, m, element_type, ap)
    
    # 1. Обчислюємо локальні матриці
    elem_stiffness = s2d.compute_element_stiffness(ELv, NLv, ap, k1=k1, k2=k2)
    elem_mass = s2d.compute_element_mass_matrix(ELv, NLv, ap) 
    # 2. Обчислити матрицю реакції
    elem_reaction = s2d.compute_element_reaction_matrix(ELv, NLv, ap, beta)
    R = s2d.assemble_global_matrix(NLv, ELv, p, m, elem_reaction, ap)
    # 2. Збираємо глобальні матриці
    K = s2d.assemble_global_matrix(NLv, ELv, p, m, elem_stiffness, ap)
    M = s2d.assemble_global_matrix(NLv, ELv, p, m, elem_mass, ap)

    # 3. Замінити K на (K + R) скрізь у схемі Кранка–Ніколсона
    KR = K + R
    
    # 3. Налаштування часу
    dt = 0.01          
    t_end = 1.0        
    num_steps = int(t_end / dt)
    theta = 0.5        # Кранк-Ніколсон

    save_every = 5   # зберігати кожен 5-й крок
    U_history = []
    time_history = []

    
    # 4. Формуємо ліву матрицю A
    A_matrix = M + theta * dt * KR
    A_matrix = f2e.apply_boundary_conditions_matrix(A_matrix, p, m, ug, ap)
    
    # 5. Початкові умови U^0 при t = 0
    # ВАЖЛИВО: Вони тепер не нульові! Ми маємо взяти точні значення функції в момент часу 0
    num_total_nodes = (ap * p + 1) * (ap * m + 1)
    U_n = np.zeros(num_total_nodes)
    for i, (x_val, y_val) in enumerate(NLv):
        U_n[i] = exact_solution_t(x_val, y_val, 0)
    
    print("Обчислення просторового вектора навантаження...")
    F_base = np.zeros(num_total_nodes)
    # Обчислюємо інтеграли Гаусса ТІЛЬКИ ОДИН РАЗ для просторової частини
    F_base = s2d.set_up_vector(fv_space_reaction, NLv, ELv, p, m, ap)
    
    print("Починаємо інтегрування по часу...")
    current_t = 0.0
    for step in range(num_steps):
        t_n = current_t
        t_n_plus_1 = current_t + dt
        
        # Масштабуємо базовий вектор на експоненту для відповідного моменту часу
        F_n = F_base * np.exp(-t_n)
        F_n_plus_1 = F_base * np.exp(-t_n_plus_1)
        
        # Права частина
        right_matrix = M - (1 - theta) * dt * KR
        b_vector = right_matrix @ U_n + dt * (theta * F_n_plus_1 + (1 - theta) * F_n)
        
        # Граничні умови
        b_vector = f2e.apply_boundary_conditions_vector(b_vector, p, m, NLv, ug, ap, current_time=t_n_plus_1)
        
        # Розв'язок
        U_next = np.linalg.solve(A_matrix, b_vector)
        U_n = U_next
        if (step + 1) % save_every == 0:
            U_history.append(U_n.copy())
            time_history.append(current_t)
        
        current_t += dt
        
        if (step + 1) % 10 == 0:
            print(f"Крок {step+1:3d}/{num_steps}, Час = {current_t:.3f}")
            
    print("Інтегрування завершено!")
    
    # 6. ВЕРИФІКАЦІЯ
    # Створюємо лямбда-функцію, яка фіксує кінцевий час t_end для передачі у твою функцію малювання
    exact_final = lambda x, y: exact_solution_t(x, y, t_end)

    print("\nАналіз похибки в кінцевий момент часу:")
    g2d.plot_2d_solution2(U_n, NLv, ELv, exact_solution=exact_final)
    g2d.animate_solution(U_history, time_history, NLv, 
                    title="Реакція-дифузія", save_filename=f"reaction_diffusion_dt_{dt}.gif")

    err_L2, seminorm_H1, err_W12 = g2d.compute_H1_error(
    exact_final, U_n, NLv, ELv, ap)
    print(f"\n--- Похибки в кінцевий момент t = {t_end} ---")
    print(f"  L2-норма:        {err_L2:.4e}")
    print(f"  H1-семінорма:    {seminorm_H1:.4e}")
    print(f"  W_2^1 (H1)-норма: {err_W12:.4e}")

    
    # test1
    # NLt1, ELt1 = m2d.uniform_mesh_with_vertices(verticest1, p, m, element_type, ap)
    # f_loadt1 = s2d.set_up_vector(ft1, NLt1, ELt1, p, m, ap)
    # elem_matricest1 = s2d.compute_element_stiffness(ELt1, NLt1, ap, k1=k1, k2=k2)
    # matrixt1 = s2d.assemble_global_stiffness_matrix(NLt1, ELt1, p, m, elem_matricest1, ap)
    # matrixt1, f_loadt1 = f2e.apply_boundary_conditions(matrixt1, f_loadt1, p, m, NLt1, ug, ap)
    # ut1 = np.linalg.solve(matrixt1, f_loadt1)
    # g2d.plot_2d_solution(ut1, NLt1, ELt1)
    # print("Координати вузлів:\n", NLt1)
    # print("Елементи:\n", ELt1)
    # print(f_loadt1)
    # print(matrixt1)
    # print(ut1)

    # test2
    # NLt2, ELt2= m2d.uniform_mesh_with_vertices(verticest2, p, m, element_type, ap)
    # f_loadt2 = s2d.set_up_vector(ft2, NLt2, ELt2, p, m, ap)
    # elem_matricest2 = s2d.compute_element_stiffness(ELt2, NLt2, ap, k1=k1, k2=k2)
    # matrixt2 =  s2d.assemble_global_stiffness_matrix(NLt2, ELt2, p, m, elem_matricest2, ap)
    # matrixt2, f_loadt2 = f2e.apply_boundary_conditions(matrixt2, f_loadt2, p, m, NLt2, ug, ap)
    # ut2 = np.linalg.solve(matrixt2, f_loadt2)
    # g2d.plot_2d_solution(ut2, NLt2, ELt2)
    # print("Координати вузлів:\n", NLt2)
    # print("Елементи:\n", ELt2)
    # print(f_loadt2)
    # print(matrixt2)
    # print(ut2)

    # test3
    # NLt3, ELt3= m2d.uniform_mesh_with_vertices(verticest3, p, m, element_type, ap)
    # f_loadt3 = s2d.set_up_vector_point_sources(sourcest3, strengthst3, NLt3, p, m, ap)
    # elem_matricest3 = s2d.compute_element_stiffness(ELt3, NLt3, ap, k1=k1, k2=k2)
    # matrixt3 =  s2d.assemble_global_stiffness_matrix(NLt3, ELt3, p, m, elem_matricest3, ap)
    # matrixt3, f_loadt3 = f2e.apply_boundary_conditions(matrixt3, f_loadt3, p, m, NLt3, ug, ap)
    # ut3 = np.linalg.solve(matrixt3, f_loadt3)
    # g2d.plot_2d_solution(ut3, NLt3, ELt3)
    # print("Координати вузлів:\n", NLt3)
    # print("Елементи:\n", ELt3)
    # print(f_loadt3)
    # print(matrixt3)
    # print(ut3)

    # test4
    # NLt4, ELt4= m2d.uniform_mesh_with_vertices(verticest4, p, m, element_type, ap)
    # f_loadt4 = s2d.set_up_vector_point_sources(sourcest4, strengthst4, NLt4, p, m, ap)
    # elem_matricest4 = s2d.compute_element_stiffness(ELt4, NLt4, ap, k1=k1, k2=k2)
    # matrixt4 =  s2d.assemble_global_stiffness_matrix(NLt4, ELt4, p, m, elem_matricest4, ap)
    # matrixt4, f_loadt4 = f2e.apply_boundary_conditions(matrixt4, f_loadt4, p, m, NLt4, ug, ap)
    # ut4 = np.linalg.solve(matrixt4, f_loadt4)
    # g2d.plot_2d_solution(ut4, NLt4, ELt4)
    # print("Координати вузлів:\n", NLt4)
    # print("Елементи:\n", ELt4)
    # print(f_loadt4)
    # print(matrixt4)
    # print(ut4)

    # test5
    # NLt5, ELt5= m2d.uniform_mesh_with_vertices(verticest5, p, m, element_type, ap)
    # f_loadt5 = s2d.set_up_vector(ft5, NLt5, ELt5, p, m, ap)
    # elem_matricest5 = s2d.compute_element_stiffness(ELt5, NLt5, ap, k1=k1t5, k2=k2t5)
    # matrixt5 =  s2d.assemble_global_stiffness_matrix(NLt5, ELt5, p, m, elem_matricest5, ap)
    # matrixt5, f_loadt5 = f2e.apply_boundary_conditions(matrixt5, f_loadt5, p, m, NLt5, ug, ap)
    # ut5 = np.linalg.solve(matrixt5, f_loadt5)
    # g2d.plot_2d_solution(ut5, NLt5, ELt5)
    # print("Координати вузлів:\n", NLt5)
    # print("Елементи:\n", ELt5)
    # print(f_loadt5)
    # print(matrixt5)
    # print(ut5)

    # test6
    # NLt6, ELt6= m2d.uniform_mesh_with_vertices(verticest6, p, m, element_type, ap)
    # f_loadt6 = s2d.set_up_vector(ft6, NLt6, ELt6, p, m, ap)
    # elem_matricest6 = s2d.compute_element_stiffness(ELt6, NLt6, ap, k1=k1t6, k2=k2t6)
    # matrixt6 =  s2d.assemble_global_stiffness_matrix(NLt6, ELt6, p, m, elem_matricest6, ap)
    # matrixt6, f_loadt6 = f2e.apply_boundary_conditions(matrixt6, f_loadt6, p, m, NLt6, ug, ap)
    # ut6 = np.linalg.solve(matrixt6, f_loadt6)
    # g2d.plot_2d_solution(ut6, NLt6, ELt6)
    # print("Координати вузлів:\n", NLt6)
    # print("Елементи:\n", ELt6)
    # print(f_loadt6)
    # print(matrixt6)
    # print(ut6)

if __name__ == '__main__':
    main()
