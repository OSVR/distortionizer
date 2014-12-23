#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDesktopWidget>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // Make the central widget have no margins, so we fill the whole region.
    QMainWindow::centralWidget()->layout()->setContentsMargins(0, 0, 0, 0);

    // Display full resolution on the second screen.
    QRect screenres = QApplication::desktop()->screenGeometry(1/*screenNumber*/);
    this->move(QPoint(screenres.x(), screenres.y()));
    this->resize(screenres.width(), screenres.height());
    this->showFullScreen();
}

MainWindow::~MainWindow()
{
    delete ui;
}
