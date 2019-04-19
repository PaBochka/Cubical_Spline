#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <vector>
#include "/home/pavel/QTproj/Spline/splyne.h"
#include <iostream>
#include <cmath>
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->page->setInteraction(QCP::iRangeZoom, true);
    ui->page->setInteraction(QCP::iRangeDrag, true);
    ui->page->axisRect()->setRangeZoom(Qt::Vertical);
    ui->page->axisRect()->setRangeDrag(Qt::Horizontal);
    ui->page->axisRect()->setRangeDrag(Qt::Vertical);
    ui->page_2->setInteraction(QCP::iRangeZoom, true);
    ui->page_2->setInteraction(QCP::iRangeDrag, true);
    ui->page_2->axisRect()->setRangeZoom(Qt::Vertical);
    ui->page_2->axisRect()->setRangeDrag(Qt::Vertical);
    ui->page_2->axisRect()->setRangeDrag(Qt::Horizontal);

    ui->page_3->setInteraction(QCP::iRangeZoom, true);
    ui->page_3->setInteraction(QCP::iRangeDrag, true);
    ui->page_3->axisRect()->setRangeZoom(Qt::Vertical);
    ui->page_3->axisRect()->setRangeDrag(Qt::Vertical);
    ui->page_3->axisRect()->setRangeDrag(Qt::Horizontal);

    connect(ui->comboBox, SIGNAL(activated(int)), ui->stackedWidget, SLOT(setCurrentIndex(int)));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    unsigned int flag = 0;
    flag = static_cast<unsigned int>(ui->comboBox_2->currentIndex());
    Splyne *splyne = new Splyne();
    int n;
    int a, b;
    int N;
    double MaxFuncError = 0;
    double xMaxFuncError = 0;
    double xMaxDev1Error = 0;
    double MaxDev1Error = 0;
    double xMaxDev2Error = 0;
    double MaxDev2Error = 0;
    QString str = ui->lineEdit->text();
    n = str.split(" ")[0].toInt();
    if (flag ==0)
    {
        a = -1;
        b = 1;
    }
    else
    {
        a = 0;
        b = 1;
    }
    double h = double(b - a) / double(n);
    double myu1, myu2;
    myu1 = ui->lineEdit_2->text().split(" ")[0].toDouble();
    myu2 = ui->lineEdit_3->text().split(" ")[0].toDouble();
//   myu1 = splyne->second_dev_fi(-1, flag);
//   myu2 = splyne->second_dev_fi(1, flag);
    std::vector<double> A(n + 1);
    std::vector<double> B(n + 1);
    std::vector<double> C(n + 1);
    std::vector<double> D(n + 1);
    std::vector<double> X(n + 1);
    for (int i = 0; i <= n; i++)
    {
        X[i] = a + i * h;
        A[i] = splyne->func_fi(X[i], flag);
    }
    C = splyne->TDMASolve(myu1, myu2, A, C, n, h, a, b);
    for (int i = 1; i <= n; i++)
    {
        B[i] = (A[i] - A[i - 1]) / h + h * (2 * C[i] + C[i - 1]) / 6;
        D[i] = (C[i] - C[i - 1]) / h;
    }

    QVector<double> y_acc;
    //y_acc = QVector<double>::fromStdVector(A);
    QVector<double> x;
    //x = QVector<double>::fromStdVector(X);
    QVector<double> y_num;
    QVector<double> error;
    QVector<double> y_dev1_acc;
    QVector<double> y_dev1_num;
    QVector<double> error_dev1;
    QVector<double> y_dev2_acc;
    QVector<double> y_dev2_num;
    QVector<double> error_dev2;
    for (double val = a; val < b; val += h)
    {
        x.push_back(val);
        y_acc.push_back(splyne->func_fi(val, flag));
        y_num.push_back(splyne->spline_s(A, D, C, B, val, X, n));
        error.push_back(splyne->func_fi(val, flag) - splyne->spline_s(A, D, C, B, val, X, n));
        y_dev1_acc.push_back(splyne->first_dev_fi(val, flag));
        y_dev2_acc.push_back(splyne->second_dev_fi(val, flag));
        y_dev1_num.push_back(splyne->first_dev_spline_s(D, C, B, val, X, n));
        y_dev2_num.push_back(splyne->second_dev_spline_s(D, C, val, X, n));
        error_dev1.push_back(splyne->first_dev_fi(val, flag) - splyne->first_dev_spline_s(D, C, B, val, X, n));
        error_dev2.push_back(splyne->second_dev_fi(val, flag) - splyne->second_dev_spline_s(D, C, val, X, n));
    }
    x.push_back(b);
    y_acc.push_back(splyne->func_fi(b, flag));
    y_num.push_back(splyne->spline_s(A, D, C, B, b, X, n));
    error.push_back(splyne->func_fi(b, flag) - splyne->spline_s(A, D, C, B, b, X, n));
    y_dev1_acc.push_back(splyne->first_dev_fi(b, flag));
    y_dev2_acc.push_back(splyne->second_dev_fi(b, flag));
    y_dev1_num.push_back(splyne->first_dev_spline_s(D, C, B, b, X, n));
    y_dev2_num.push_back(splyne->second_dev_spline_s(D, C, b, X, n));
    error_dev1.push_back(splyne->first_dev_fi(b, flag) - splyne->first_dev_spline_s(D, C, B, b, X, n));
    error_dev2.push_back(splyne->second_dev_fi(b, flag) - splyne->second_dev_spline_s(D, C, b, X, n));
    //График функции и сплайна
    ui->page->clearGraphs();
    ui->page->addGraph();
    ui->page->addGraph();
    ui->page->addGraph();
    ui->page->graph(0)->setData(x, y_acc);
    ui->page->graph(0)->setPen(QPen(Qt::red));
    ui->page->graph(1)->setData(x, y_num);
    ui->page->graph(1)->setPen(QPen(Qt::blue));
    ui->page->graph(2)->setData(x, error);
    ui->page->graph(2)->setPen(QPen(Qt::green));
    ui->page->rescaleAxes();
    ui->page->xAxis->setLabel("x");
    ui->page->yAxis->setLabel("y");
    ui->page->legend->setVisible(true);
    ui->page->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft|Qt::AlignTop);
    ui->page->graph(0)->setName("Accuracy");
    ui->page->graph(1)->setName("Number");
    ui->page->graph(2)->setName("accuracy - number");
    ui->page->graph(0)->setAntialiased(false);
    ui->page->replot();


    //График первой производной
    ui->page_2->clearGraphs();
    ui->page_2->addGraph();
    ui->page_2->addGraph();
    ui->page_2->addGraph();
    ui->page_2->graph(0)->setData(x, y_dev1_acc);
    ui->page_2->graph(0)->setPen(QPen(Qt::red));
    ui->page_2->graph(1)->setData(x, y_dev1_num);
    ui->page_2->graph(1)->setPen(QPen(Qt::blue));
    ui->page_2->graph(2)->setData(x, error_dev1);
    ui->page_2->graph(2)->setPen(QPen(Qt::green));
    ui->page_2->rescaleAxes();
    ui->page_2->xAxis->setLabel("x");
    ui->page_2->yAxis->setLabel("y");
    ui->page_2->legend->setVisible(true);
    ui->page_2->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft|Qt::AlignTop);
    ui->page_2->graph(0)->setName("Accuracy");
    ui->page_2->graph(1)->setName("Number");
    ui->page_2->graph(2)->setName("accuracy - number");
    ui->page_2->graph(0)->setAntialiased(false);
    ui->page_2->replot();


    //График второй производной
    ui->page_3->clearGraphs();
    ui->page_3->addGraph();
    ui->page_3->addGraph();
    ui->page_3->addGraph();
    ui->page_3->graph(0)->setData(x, y_dev2_acc);
    ui->page_3->graph(0)->setPen(QPen(Qt::red));
    ui->page_3->graph(1)->setData(x, y_dev2_num);
    ui->page_3->graph(1)->setPen(QPen(Qt::blue));
    ui->page_3->graph(2)->setData(x, error_dev2);
    ui->page_3->graph(2)->setPen(QPen(Qt::green));
    ui->page_3->rescaleAxes();
    ui->page_3->xAxis->setLabel("x");
    ui->page_3->yAxis->setLabel("y");
    ui->page_3->legend->setVisible(true);
    ui->page_3->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft|Qt::AlignTop);
    ui->page_3->graph(0)->setName("Accuracy");
    ui->page_3->graph(1)->setName("Number");
    ui->page_3->graph(2)->setName("accuracy - number");
    ui->page_3->graph(0)->setAntialiased(false);
    ui->page_3->replot();

    for (double val = a; val < b; val += h / 4.0)
    {
        auto tmp = splyne->func_fi(val, flag) - splyne->spline_s(A, D, C, B, val, X, n);
        auto tmp_dev1 = splyne->first_dev_fi(val,flag) - splyne->first_dev_spline_s(D, C, B, val, X, n);
        auto tmp_dev2 = splyne->second_dev_fi(val, flag) - splyne->second_dev_spline_s(D, C, val, X, n);
        if(MaxFuncError < tmp)
        {
            MaxFuncError = tmp;
            xMaxFuncError = val;
        }
        if(MaxDev1Error < tmp_dev1)
        {
            MaxDev1Error = tmp_dev1;
            xMaxDev1Error = val;
        }
        if(MaxDev2Error < tmp_dev2)
        {
            MaxDev2Error = tmp_dev2;
            xMaxDev2Error = val;
        }

    }
    if(n < 100)
        N = n;
    else
        N = 100;
    ui->tableWidget->setRowCount(N + 1);
    for (int i = 1; i <= N; i++) {
        ui->tableWidget->setItem(i, 0, new QTableWidgetItem(QString::number(i)));
        ui->tableWidget->setItem(i, 1, new QTableWidgetItem(QString::number(X[i - 1], 'f', 7)));
        ui->tableWidget->setItem(i, 2, new QTableWidgetItem(QString::number(X[i], 'f', 7)));
        ui->tableWidget->setItem(i, 3, new QTableWidgetItem(QString::number(A[i], 'f', 7)));
        ui->tableWidget->setItem(i, 4, new QTableWidgetItem(QString::number(B[i], 'f', 7)));
        ui->tableWidget->setItem(i, 5, new QTableWidgetItem(QString::number(C[i], 'f', 7)));
        ui->tableWidget->setItem(i, 6, new QTableWidgetItem(QString::number(D[i], 'f', 7)));
    }
    if(n < 100)
        N = 4 * n;
    else
        N = 100;
    ui->tableWidget_2->setRowCount(N + 1);
    double y = a;
    for (int i = 0; i <= N; i++) {
        ui->tableWidget_2->setItem(i, 0, new QTableWidgetItem(QString::number(i)));
        ui->tableWidget_2->setItem(i, 1, new QTableWidgetItem(QString::number(y, 'f', 7)));
        ui->tableWidget_2->setItem(i, 2, new QTableWidgetItem(QString::number(splyne->func_fi(y, flag), 'f', 7)));
        ui->tableWidget_2->setItem(i, 3, new QTableWidgetItem(QString::number(splyne->spline_s(A, D, C, B, y, X, n), 'f', 7)));
        ui->tableWidget_2->setItem(i, 4, new QTableWidgetItem(QString::number(abs(splyne->func_fi(y, flag) - splyne->spline_s(A, D, C, B, y, X, n)), 'f', 7)));
        ui->tableWidget_2->setItem(i, 5, new QTableWidgetItem(QString::number(splyne->first_dev_fi(y, flag), 'f', 7)));
        ui->tableWidget_2->setItem(i, 6, new QTableWidgetItem(QString::number(splyne->first_dev_spline_s(D, C, B, y, X, n), 'f', 7)));
        ui->tableWidget_2->setItem(i, 7, new QTableWidgetItem(QString::number(abs(splyne->first_dev_fi(y, flag) - splyne->first_dev_spline_s(D, C, B, y, X, n)), 'f', 7)));
        y += h / 4;
        y = round(y*100000)/100000;

    }
    ui->tableWidget_3->setRowCount(N + 1);
    y = a;
    for (int i = 0; i <= N; i++) {
        ui->tableWidget_3->setItem(i, 0, new QTableWidgetItem(QString::number(i)));
        ui->tableWidget_3->setItem(i, 1, new QTableWidgetItem(QString::number(y, 'f')));
        ui->tableWidget_3->setItem(i, 2, new QTableWidgetItem(QString::number(splyne->second_dev_fi(y, flag), 'f')));
        ui->tableWidget_3->setItem(i, 3, new QTableWidgetItem(QString::number(splyne->second_dev_spline_s(D, C, y, X, n), 'f')));
        ui->tableWidget_3->setItem(i, 4, new QTableWidgetItem(QString::number(abs(splyne->second_dev_fi(y, flag) - splyne->second_dev_spline_s(D, C, y, X, n)), 'f')));
        y += h / 4;
        y = round(y*100000)/100000;
    }
    ui->textBrowser->setText(" Cетка сплайна: n = " + QString::number(n) + " \n Контрольная сетка N = " + QString::number(N)
                             + " \n max|F(x) - S(x)| = " + QString::number(MaxFuncError, 'f', 16) + " при x = " +
                             QString::number(xMaxFuncError, 'f', 4) + " \n max|F'(x) - S'(x)| = "
                             + QString::number(MaxDev1Error, 'f', 16) + " при x = " +
                             QString::number(xMaxDev1Error, 'f', 4) + " \n max|F''(x) - S''(x)| = "
                             + QString::number(MaxDev2Error, 'f', 16) + " при x = " +
                             QString::number(xMaxDev2Error, 'f', 4)+ ".");
    delete splyne;

};



void MainWindow::on_pushButton_2_clicked()
{
    Splyne *splyne = new Splyne();
    unsigned int flag = static_cast<unsigned int>(ui->comboBox_2->currentIndex());
    int a, b;
    if(flag == 0)
    {
        a = -1;
        b = 1;
    }
    else
    {
        a = 0;
        b = 1;
    }

    auto myu1 = splyne->second_dev_fi(a, flag);
    auto myu2 = splyne->second_dev_fi(b, flag);

    ui->lineEdit_2->setText(QString::number(myu1));
    ui->lineEdit_3->setText(QString::number(myu2));
    delete splyne;
}
