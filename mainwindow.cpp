#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <vector>
#include "/home/pavel/QTproj/Spline/splyne.h"
#include <iostream>
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->widget->setInteraction(QCP::iRangeZoom, true);
    ui->widget->setInteraction(QCP::iRangeDrag, true);
    ui->widget->axisRect()->setRangeZoom(Qt::Vertical);
    ui->widget->axisRect()->setRangeDrag(Qt::Vertical);

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    Splyne *splyne = new Splyne();
    int n, a, b;
    QString str = ui->lineEdit->text();
    n = str.split(" ")[0].toInt();
    a = -1;
    b = 1;
    double h = double(b - a) / double(n);
    double myu1, myu2;
    myu1 = splyne->second_dev_fi(-1);
    myu2 = splyne->second_dev_fi(1);
    std::vector<double> A(n + 1);
    std::vector<double> B(n + 1);
    std::vector<double> C(n + 1);
    std::vector<double> D(n + 1);
    std::vector<double> X(n + 1);
    for (int i = 0; i <= n; i++)
    {
        X[i] = a + i * h;
        A[i] = splyne->func_fi(X[i]);
    }
    C = splyne->TDMASolve(myu1, myu2, A, C, n, h, a, b);
    for (int i = 1; i <= n; i++)
    {
        B[i] = (A[i] - A[i - 1]) / h + h * (2 * C[i] + C[i - 1]) / 6;
        D[i] = (C[i] - C[i - 1]) / h;
    }

    QVector<double> y_acc(n + 1);
    //y_acc = QVector<double>::fromStdVector(A);
    QVector<double> x(n + 1);
    //x = QVector<double>::fromStdVector(X);
    QVector<double> y_num(n + 1);
    int i = 0;
    for (double val = a; (val < b || i <= n); val += h)
    {
        x.push_back(val);
        y_acc.push_back(splyne->func_fi(val));
        y_num.push_back(splyne->spline_s(A[i], D[i], C[i], B[i], val, X[i]));
        i++;
    }
    ui->widget->clearGraphs();
    ui->widget->addGraph();
    ui->widget->addGraph();
    ui->widget->graph(0)->setData(x, y_acc);
    ui->widget->graph(0)->setPen(QPen(Qt::red));
    ui->widget->graph(1)->setData(x, y_num);
    ui->widget->graph(1)->setPen(QPen(Qt::blue));
    ui->widget->rescaleAxes();
    ui->widget->xAxis->setLabel("x");
    ui->widget->yAxis->setLabel("y");
    ui->widget->legend->setVisible(true);
    ui->widget->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft|Qt::AlignTop);
    ui->widget->graph(0)->setName("Accuracy");
    ui->widget->graph(1)->setName("Number");
    ui->widget->graph(0)->setAntialiased(false);


    ui->widget->replot();
}


