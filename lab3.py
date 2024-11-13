import numpy as numpy
import matplotlib.pyplot as matplotlib


km = 10^-7

size = 81
zero_index = int((size - 1) / 2)

axis = numpy.linspace(- (size - 1) / 2, (size - 1) / 2, size)
X, Y, Z = numpy.meshgrid(axis, axis, axis)

I = numpy.zeros(X.shape)

Dx = numpy.zeros(X.shape)
Dy = numpy.zeros(Y.shape)
Dz = numpy.zeros(Z.shape)

Bx = numpy.zeros(X.shape)
By = numpy.zeros(Y.shape)
Bz = numpy.zeros(Z.shape)


# Recta
def generar_recta(corriente, largo):
    k1 = int(zero_index - largo / 2)
    k2 = int(zero_index + largo / 2)
    for k in range(k1, k2):
        Dz[zero_index][zero_index][k] = 1
        I[zero_index][zero_index][k] = corriente

# Espira
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


""" generar_recta(10, 81) """
""" generar_espira(50, 15, 15) """
""" generar_solenoide(200, 10, 15, 41, 8) """
""" generar_bobina_helmholtz(50, 20, 15) """


Dl = numpy.sqrt(numpy.square(Dx) + numpy.square(Dy) + numpy.square(Dz))
conductores = Dl != 0

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

    r[r == 0] = numpy.inf

    k = km * corriente / r**3

    Bx += k * (dy * dist_z - dz * dist_y)
    By += k * (dz * dist_x - dx * dist_z)
    Bz += k * (dx * dist_y - dy * dist_x)

    
calculate_load = input("Enter 0 to calculate | 1 to load last field: ")


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


chosen_plane = input("Select a plane: xy: 0 | yz: 1 | zx: 2: ")
third_axis_value = input("Select the position of the chosen plane in the third axis: ")

if chosen_plane == "0":
    first_axis = X[:,:,int(third_axis_value) + zero_index]
    second_axis = Y[:,:,int(third_axis_value) + zero_index]
    first_axis_field_components = Bx[:,:,int(third_axis_value) + zero_index]
    second_axis_field_components = By[:,:,int(third_axis_value) + zero_index]
    conductors_mask = conductores[:,:,int(third_axis_value) + zero_index]
elif chosen_plane == "1":
    first_axis = numpy.transpose(Y[:,int(third_axis_value) + zero_index,:])
    second_axis = numpy.transpose(Z[:,int(third_axis_value) + zero_index,:])
    first_axis_field_components = numpy.transpose(By[:,int(third_axis_value) + zero_index,:])
    second_axis_field_components = numpy.transpose(Bz[:,int(third_axis_value) + zero_index,:])
    conductors_mask = numpy.transpose(conductores[:,int(third_axis_value) + zero_index,:])
elif chosen_plane == "2":
    first_axis = Z[int(third_axis_value) + zero_index,:,:]
    second_axis = X[int(third_axis_value) + zero_index,:,:]
    first_axis_field_components = Bz[int(third_axis_value) + zero_index,:,:]
    second_axis_field_components = Bx[int(third_axis_value) + zero_index,:,:]
    conductors_mask = conductores[int(third_axis_value) + zero_index,:,:]
else:
    quit()


plane_field_magnitude = numpy.sqrt(numpy.square(first_axis_field_components) + numpy.square(second_axis_field_components))


fig = matplotlib.figure(figsize=matplotlib.figaspect(0.4))


ax = fig.add_subplot(1, 2, 1)

ax.scatter(first_axis[conductors_mask], second_axis[conductors_mask], color='dimgray', s=100)
ax.streamplot(first_axis, second_axis, first_axis_field_components, second_axis_field_components, color='cornflowerblue', linewidth=1, density=2)

ax.set_xlim(- (size - 1) / 2, (size - 1) / 2)
ax.set_ylim(- (size - 1) / 2, (size - 1) / 2)


ax3d = fig.add_subplot(1, 2, 2, projection='3d')
surf = ax3d.plot_surface(first_axis, second_axis, plane_field_magnitude, cmap='Blues', linewidth=0, antialiased=True)


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
fig.colorbar(surf, ax=ax3d, orientation='vertical', label='Magnetic Field Magnitude [T]')

# Datos extra para el grafico 2
ax3d.set_title('Magnetic Field Surface')
if chosen_plane == "0":
    ax3d.set_xlabel('x [m]')
    ax3d.set_ylabel('y [m]')
    ax3d.set_zlabel('|Bxy| [T]')
if chosen_plane == "1":
    ax3d.set_xlabel('y [m]')
    ax3d.set_ylabel('z [m]')
    ax3d.set_zlabel('|Byz| [T]')
if chosen_plane == "2":
    ax3d.set_xlabel('z [m]')
    ax3d.set_ylabel('x [m]')
    ax3d.set_zlabel('|Bzx| [T]')


matplotlib.tight_layout()
matplotlib.show()