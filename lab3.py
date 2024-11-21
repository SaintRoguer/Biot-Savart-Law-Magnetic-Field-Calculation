import numpy as numpy
import matplotlib.pyplot as matplotlib

# constante igual a mu 0 sobre 4 pi
km = 10**(-7)

# rango para cada eje de coordenadas
size = 81
# indice de la matriz que se corresponde con el origen de cada eje de coordenadas
zero_index = int((size - 1) / 2)

# definicion de la region del espacio
axis = numpy.linspace(- (size - 1) / 2, (size - 1) / 2, size)
X, Y, Z = numpy.meshgrid(axis, axis, axis)

# valor de la corriente en cada punto del espacio, inicializado en ceros
I = numpy.zeros(X.shape)

# valores de las componentes del vector diferencial de linea en cada punto del espacio, inicializado en ceros
Dx = numpy.zeros(X.shape)
Dy = numpy.zeros(Y.shape)
Dz = numpy.zeros(Z.shape)

# valores de las componentes del vector campo magnetico en cada punto del espacio, inicializado en ceros
Bx = numpy.zeros(X.shape)
By = numpy.zeros(Y.shape)
Bz = numpy.zeros(Z.shape)

# ***********************************************************************************************************
# ************************************** GENERACION DE CONDUCTORES ******************************************
# ***********************************************************************************************************

# funcion que genera un conductor recto
def generar_recta(corriente, largo):
    k1 = int(zero_index - largo / 2)
    k2 = int(zero_index + largo / 2)
    for k in range(k1, k2):
        Dz[zero_index][zero_index][k] = 1
        I[zero_index][zero_index][k] = corriente

# funcion que genera un conductor con forma de espira circular
def generar_espira(nro_diferenciales, corriente, radio):
    for t in range(0, nro_diferenciales):
        theta = t * 2 * numpy.pi / nro_diferenciales

        x = radio * numpy.cos(theta)
        y = radio * numpy.sin(theta)

        magnitud_diferencial = 2 * numpy.pi * radio / nro_diferenciales

        i = round(x) + zero_index
        j = round(y) + zero_index

        Dx[j][i][zero_index] += -y / radio * magnitud_diferencial
        Dy[j][i][zero_index] += x / radio * magnitud_diferencial
        I[j][i][zero_index] = corriente

# funcion que genera un solenoide
def generar_solenoide(nro_diferenciales, corriente, radio, largo, vueltas):
    for t in range(0, nro_diferenciales):
        theta = t * 2 * numpy.pi * vueltas / nro_diferenciales

        x = radio * numpy.cos(theta)
        y = radio * numpy.sin(theta)
        z = t * largo / nro_diferenciales - largo / 2

        longitud_solenoide = numpy.sqrt((2 * numpy.pi * radio * vueltas)**2 + largo**2)
        magnitud_diferencial = longitud_solenoide / nro_diferenciales

        i = round(x) + zero_index
        j = round(y) + zero_index
        k = round(z) + zero_index

        dxdt = -y
        dydt = x
        dzdt = largo / nro_diferenciales
        modulo_vector_tangente = numpy.sqrt(radio**2 + dzdt**2)

        Dx[j][i][k] += dxdt / modulo_vector_tangente * magnitud_diferencial
        Dy[j][i][k] += dydt / modulo_vector_tangente * magnitud_diferencial
        Dz[j][i][k] += dzdt / modulo_vector_tangente * magnitud_diferencial
        I[j][i][k] = corriente

# funcion que genera una bobina de Helmholtz
def generar_bobina_helmholtz(nro_diferenciales, corriente, radio):
    z1 = int(radio / 2)
    z2 = int(-radio / 2)
    k1 = round(z1) + zero_index
    k2 = round(z2) + zero_index
    for t in range(0, nro_diferenciales):
        theta = t * 2 * numpy.pi / nro_diferenciales

        x = radio * numpy.cos(theta)
        y = radio * numpy.sin(theta)

        magnitud_diferencial = 2 * numpy.pi * radio / nro_diferenciales

        i = round(x) + zero_index
        j = round(y) + zero_index
        
        Dx[j][i][k1] += -y / radio * magnitud_diferencial
        Dy[j][i][k1] += x / radio * magnitud_diferencial
        I[j][i][k1] = corriente

        Dx[j][i][k2] += -y / radio * magnitud_diferencial
        Dy[j][i][k2] += x / radio * magnitud_diferencial
        I[j][i][k2] = corriente


# ***********************************************************************************************************
# ************************************** PROGRAMA PRINCIPAL *************************************************
# ***********************************************************************************************************

# llamadas a funciones para la generacion de la configuracion de conductores deseada
generar_recta(10, 81)
""" generar_espira(50, 15, 15) """
""" generar_solenoide(200, 10, 15, 41, 8) """
""" generar_bobina_helmholtz(50, 20, 15) """

# Obtencion de los puntos donde el diferencial de linea es no nulo, es decir, los puntos correspondientes a los conductores
Dl = numpy.sqrt(numpy.square(Dx) + numpy.square(Dy) + numpy.square(Dz))
conductores = Dl != 0

# funcion que calcula la contribucion al campo magnetico por el elemento de corriente en un punto
def add_magnetic_field_diferential(i, j, k):
    global Bx, By, Bz

    corriente = I[j][i][k]
    dx = Dx[j][i][k]
    dy = Dy[j][i][k]
    dz = Dz[j][i][k]

    dist_x = X - (i - zero_index)
    dist_y = Y - (j - zero_index)
    dist_z = Z - (k - zero_index)
    r = numpy.sqrt(dist_x**2 + dist_y**2 + dist_z**2)

    with numpy.errstate(divide='ignore', invalid='ignore'):
        k = km * corriente / r**3

        Bx += k * (dy * dist_z - dz * dist_y)
        By += k * (dz * dist_x - dx * dist_z)
        Bz += k * (dx * dist_y - dy * dist_x)


# variable para seleccionar si se quiere calcular el campo, o usar el campo guardado en la ultima ejecucion del programa
calculate_load = input("Enter 0 to calculate | 1 to load last field: ")


# calcula el campo recorriendo todo punto de los conductores y computando su contribucion al campo y luego lo salva, o bien carga el campo guardado
if calculate_load == "0":
    for i in range(size):
        for j in range(size):
            for k in range(size):
                if conductores[j][i][k]:
                    add_magnetic_field_diferential(i, j, k)
    numpy.save("saved_fields/saved_field_x.npy", Bx)
    numpy.save("saved_fields/saved_field_y.npy", By)
    numpy.save("saved_fields/saved_field_z.npy", Bz)
elif calculate_load == "1":
    Bx = numpy.load("saved_fields/saved_field_x.npy")
    By = numpy.load("saved_fields/saved_field_y.npy")
    Bz = numpy.load("saved_fields/saved_field_z.npy")
else:
    quit()


# ***********************************************************************************************************
# *************************************** INPUT *************************************************************
# ***********************************************************************************************************

# input para el plano sobre el cual graficar el campo en 2D
chosen_plane = input("Select a plane: xy: 0 | yz: 1 | zx: 2: ")
third_axis_value = int(input("Select the position of the chosen plane in the third axis: "))

# input para la region del espacio donde graficar en 3D
print("Select the center point for the 3D graphic: ")
graphic_3d_center_x = int(input("x: "))
graphic_3d_center_y = int(input("y: "))
graphic_3d_center_z = int(input("z: "))

# input para el punto donde obtener la expresion del vector campo magnetico
print("Select a point to calculate the magnetic field vector: ")
point_x = int(input("x: "))
point_y = int(input("y: "))
point_z = int(input("z: "))

# ***********************************************************************************************************
# ************************************** OUTPUT *************************************************************
# ***********************************************************************************************************

# obtiene los valores para el punto especificado
x_value = Bx[zero_index + point_y][zero_index + point_x][zero_index + point_z]
y_value = By[zero_index + point_y][zero_index + point_x][zero_index + point_z]
z_value = Bz[zero_index + point_y][zero_index + point_x][zero_index + point_z]

# imprime el vector campo magnetico en el punto especificado
print(f"Field in ({point_x}, {point_y}, {point_z}) is: ({x_value})i + ({y_value})j + ({z_value})k")

# asignacion de variables correspondientes segun el plano elegido para la posterior grafica
if chosen_plane == "0":
    first_axis = X[:,:,third_axis_value + zero_index]
    second_axis = Y[:,:,third_axis_value + zero_index]
    first_axis_field_components = Bx[:,:,third_axis_value + zero_index]
    second_axis_field_components = By[:,:,third_axis_value + zero_index]
    conductors_mask = conductores[:,:,third_axis_value + zero_index]
elif chosen_plane == "1":
    first_axis = numpy.transpose(Y[:,third_axis_value + zero_index,:])
    second_axis = numpy.transpose(Z[:,third_axis_value + zero_index,:])
    first_axis_field_components = numpy.transpose(By[:,third_axis_value + zero_index,:])
    second_axis_field_components = numpy.transpose(Bz[:,third_axis_value + zero_index,:])
    conductors_mask = numpy.transpose(conductores[:,third_axis_value + zero_index,:])
elif chosen_plane == "2":
    first_axis = Z[third_axis_value + zero_index,:,:]
    second_axis = X[third_axis_value + zero_index,:,:]
    first_axis_field_components = Bz[third_axis_value + zero_index,:,:]
    second_axis_field_components = Bx[third_axis_value + zero_index,:,:]
    conductors_mask = conductores[third_axis_value + zero_index,:,:]
else:
    quit()

# obtencion de la magnitud del campo paralela al plano elegido
plane_field_magnitude = numpy.sqrt(numpy.square(first_axis_field_components) + numpy.square(second_axis_field_components))

# definicion de la figura para los graficos del plano
fig = matplotlib.figure(figsize=matplotlib.figaspect(0.4))
ax = fig.add_subplot(1, 2, 1)
ax_surface = fig.add_subplot(1, 2, 2, projection='3d')

# definicion de la figura para el grafico de la region del espacio
fig3d, ax3d = matplotlib.subplots(subplot_kw={"projection": "3d"})

# el primer grafico contiene un grafico de flujo que indica la direccion y sentido del campo en el plano, y un grafico de dispersion que indica los puntos correspondientes a conductores
ax.scatter(first_axis[conductors_mask], second_axis[conductors_mask], color='dimgray', s=100)
ax.streamplot(first_axis, second_axis, first_axis_field_components, second_axis_field_components, color='cornflowerblue', linewidth=1, density=2)

# el segundo grafico es una superficie que indica la magnitud del campo magnetico sobre el plano
surf = ax_surface.plot_surface(first_axis, second_axis, plane_field_magnitude, cmap='Blues', linewidth=0, antialiased=True)

# el tercer grafico muestra las direcciones y sentidos del campo magnetico, representado como vectores unitarios en cada punto de la region considerada
graphic_3d_mask = numpy.zeros_like(X, dtype=bool)
graphic_3d_mask[zero_index + graphic_3d_center_y - 5 : zero_index + graphic_3d_center_y + 5,
                zero_index + graphic_3d_center_x - 5 : zero_index + graphic_3d_center_x + 5,
                zero_index + graphic_3d_center_z - 5 : zero_index + graphic_3d_center_z + 5] = True
ax3d.quiver(X[graphic_3d_mask], Y[graphic_3d_mask], Z[graphic_3d_mask], Bx[graphic_3d_mask], By[graphic_3d_mask], Bz[graphic_3d_mask], normalize=True, length=0.5)

# Datos extra para el grafico 1
ax.set_title('Magnetic Field Lines')
if chosen_plane == "0":
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
elif chosen_plane == "1":
    ax.set_xlabel('y [m]')
    ax.set_ylabel('z [m]')
elif chosen_plane == "2":
    ax.set_xlabel('z [m]')
    ax.set_ylabel('x [m]')
fig.colorbar(surf, ax=ax_surface, orientation='vertical', label='Magnetic Field Magnitude [T]', pad=0.1)
ax.set_xlim(- (size - 1) / 2, (size - 1) / 2)
ax.set_ylim(- (size - 1) / 2, (size - 1) / 2)

# Datos extra para el grafico 2
ax_surface.set_title('Magnetic Field Surface')
if chosen_plane == "0":
    ax_surface.set_xlabel('x [m]')
    ax_surface.set_ylabel('y [m]')
    ax_surface.set_zlabel('|Bxy| [T]')
if chosen_plane == "1":
    ax_surface.set_xlabel('y [m]')
    ax_surface.set_ylabel('z [m]')
    ax_surface.set_zlabel('|Byz| [T]')
if chosen_plane == "2":
    ax_surface.set_xlabel('z [m]')
    ax_surface.set_ylabel('x [m]')
    ax_surface.set_zlabel('|Bzx| [T]')

# Datos extra para el grafico 3
ax3d.set_xlabel('x [m]')
ax3d.set_ylabel('y [m]')
ax3d.set_zlabel('z [m]')

# muestra los graficos
matplotlib.tight_layout()
matplotlib.show()