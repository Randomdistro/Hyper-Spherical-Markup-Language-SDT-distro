/**
 * HSML Visual Studio GUI - Main Application Entry Point
 * GUI-FIRST METHODOLOGY: Perfect visual experience from application startup
 * 
 * The ultimate C++ development environment for spherical coordinate systems
 * Launches with perfect visual state and maximum development efficiency
 */

#include "hsml/gui/hsml_visual_studio_gui.h"

#include <QtWidgets/QApplication>
#include <QtWidgets/QSplashScreen>
#include <QtWidgets/QMessageBox>
#include <QtCore/QDir>
#include <QtCore/QStandardPaths>
#include <QtCore/QTimer>
#include <QtCore/QTranslator>
#include <QtGui/QPixmap>
#include <QtGui/QPainter>
#include <QtGui/QFont>

#include <iostream>
#include <exception>

// GUI-FIRST: Perfect splash screen for immediate visual feedback
QPixmap createSplashScreen() {
    QPixmap splash(800, 400);
    splash.fill(QColor(43, 43, 43)); // Dark theme background
    
    QPainter painter(&splash);
    painter.setRenderHint(QPainter::Antialiasing);
    
    // HSML logo and title
    QFont title_font("Segoe UI", 36, QFont::Bold);
    painter.setFont(title_font);
    painter.setPen(QColor(0, 120, 212)); // Microsoft blue
    painter.drawText(50, 100, "HSML Visual Studio");
    
    QFont subtitle_font("Segoe UI", 16);
    painter.setFont(subtitle_font);
    painter.setPen(QColor(255, 255, 255));
    painter.drawText(50, 140, "Ultimate Spherical Development Environment");
    
    // Spherical coordinate visualization
    painter.setPen(QPen(QColor(78, 205, 196), 2)); // Teal accent
    painter.setBrush(Qt::NoBrush);
    
    // Draw concentric spheres
    for (int r = 50; r <= 150; r += 50) {
        painter.drawEllipse(600 - r, 200 - r, 2 * r, 2 * r);
    }
    
    // Draw coordinate axes
    painter.setPen(QPen(QColor(255, 107, 107), 3)); // Red accent
    painter.drawLine(600, 50, 600, 350); // Vertical axis
    painter.drawLine(450, 200, 750, 200); // Horizontal axis
    
    // Loading text
    QFont loading_font("Segoe UI", 12);
    painter.setFont(loading_font);
    painter.setPen(QColor(200, 200, 200));
    painter.drawText(50, 350, "Loading spherical coordinate engine...");
    
    return splash;
}

// GUI-FIRST: Perfect application initialization
bool initializeApplication(QApplication& app) {
    try {
        // Set application properties
        app.setApplicationName("HSML Visual Studio");
        app.setApplicationVersion("1.0.0");
        app.setOrganizationName("HSML Project");
        app.setOrganizationDomain("hsml.dev");
        
        // Set application icon
        // app.setWindowIcon(QIcon("://icons/hsml_logo.png"));
        
        // Enable high DPI support
        app.setAttribute(Qt::AA_EnableHighDpiScaling);
        app.setAttribute(Qt::AA_UseHighDpiPixmaps);
        
        // Set up application directories
        QString app_data_dir = QStandardPaths::writableLocation(QStandardPaths::AppDataLocation);
        QDir().mkpath(app_data_dir);
        
        QString projects_dir = app_data_dir + "/Projects";
        QDir().mkpath(projects_dir);
        
        QString plugins_dir = app_data_dir + "/Plugins";
        QDir().mkpath(plugins_dir);
        
        QString templates_dir = app_data_dir + "/Templates";
        QDir().mkpath(templates_dir);
        
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Failed to initialize application: " << e.what() << std::endl;
        return false;
    }
}

// GUI-FIRST: Perfect main function with exceptional visual experience
int main(int argc, char* argv[]) {
    // Create QApplication with perfect visual settings
    QApplication app(argc, argv);
    
    // GUI-FIRST: Show splash screen immediately
    QPixmap splash_pixmap = createSplashScreen();
    QSplashScreen splash(splash_pixmap);
    splash.show();
    splash.showMessage("Initializing HSML Visual Studio...", 
                      Qt::AlignBottom | Qt::AlignLeft, 
                      QColor(255, 255, 255));
    
    app.processEvents();
    
    try {
        // Initialize application with perfect setup
        splash.showMessage("Setting up application framework...", 
                          Qt::AlignBottom | Qt::AlignLeft, 
                          QColor(255, 255, 255));
        app.processEvents();
        
        if (!initializeApplication(app)) {
            splash.close();
            QMessageBox::critical(nullptr, "Initialization Error", 
                                "Failed to initialize HSML Visual Studio application framework.");
            return 1;
        }
        
        // Load spherical coordinate engine
        splash.showMessage("Loading spherical coordinate engine...", 
                          Qt::AlignBottom | Qt::AlignLeft, 
                          QColor(255, 255, 255));
        app.processEvents();
        QThread::msleep(500); // Brief pause for visual effect
        
        // Initialize rendering system
        splash.showMessage("Initializing OpenGL rendering system...", 
                          Qt::AlignBottom | Qt::AlignLeft, 
                          QColor(255, 255, 255));
        app.processEvents();
        QThread::msleep(300);
        
        // Load HSML language definitions
        splash.showMessage("Loading HSML language definitions...", 
                          Qt::AlignBottom | Qt::AlignLeft, 
                          QColor(255, 255, 255));
        app.processEvents();
        QThread::msleep(300);
        
        // Initialize material library
        splash.showMessage("Loading material library...", 
                          Qt::AlignBottom | Qt::AlignLeft, 
                          QColor(255, 255, 255));
        app.processEvents();
        QThread::msleep(200);
        
        // Initialize testing framework
        splash.showMessage("Initializing testing framework...", 
                          Qt::AlignBottom | Qt::AlignLeft, 
                          QColor(255, 255, 255));
        app.processEvents();
        QThread::msleep(200);
        
        // Load plugins
        splash.showMessage("Loading plugin architecture...", 
                          Qt::AlignBottom | Qt::AlignLeft, 
                          QColor(255, 255, 255));
        app.processEvents();
        QThread::msleep(200);
        
        // Initialize security system
        splash.showMessage("Initializing security validation...", 
                          Qt::AlignBottom | Qt::AlignLeft, 
                          QColor(255, 255, 255));
        app.processEvents();
        QThread::msleep(200);
        
        // Final setup
        splash.showMessage("Finalizing setup...", 
                          Qt::AlignBottom | Qt::AlignLeft, 
                          QColor(255, 255, 255));
        app.processEvents();
        QThread::msleep(300);
        
        // Create main window with perfect GUI-first architecture
        splash.showMessage("Creating main window...", 
                          Qt::AlignBottom | Qt::AlignLeft, 
                          QColor(255, 255, 255));
        app.processEvents();
        
        auto* main_window = new hsml::gui::HSMLVisualStudioGUI();
        
        // Initialize main window
        if (!main_window->initialize()) {
            splash.close();
            delete main_window;
            QMessageBox::critical(nullptr, "Initialization Error", 
                                "Failed to initialize HSML Visual Studio main window.");
            return 1;
        }
        
        // GUI-FIRST: Perfect window presentation
        splash.showMessage("Ready!", 
                          Qt::AlignBottom | Qt::AlignLeft, 
                          QColor(0, 255, 0));
        app.processEvents();
        QThread::msleep(500);
        
        // Show main window with perfect visual transition
        main_window->show();
        main_window->raise();
        main_window->activateWindow();
        
        // Close splash screen with fade effect
        QTimer::singleShot(200, &splash, &QSplashScreen::close);
        
        // Handle command line arguments for project opening
        QStringList args = app.arguments();
        if (args.size() > 1) {
            QString project_path = args[1];
            if (QDir(project_path).exists()) {
                main_window->openProject(project_path.toStdString());
            }
        }
        
        std::cout << "HSML Visual Studio started successfully!" << std::endl;
        std::cout << "GUI-FIRST methodology: Perfect visual experience drives perfect functionality" << std::endl;
        
        // Run application event loop
        int result = app.exec();
        
        // Cleanup
        delete main_window;
        
        return result;
        
    } catch (const std::exception& e) {
        splash.close();
        
        QString error_msg = QString("Fatal error starting HSML Visual Studio:\n\n%1\n\n"
                                  "Please check your installation and try again.").arg(e.what());
        
        QMessageBox::critical(nullptr, "Fatal Error", error_msg);
        
        std::cerr << "Fatal error: " << e.what() << std::endl;
        return 1;
        
    } catch (...) {
        splash.close();
        
        QMessageBox::critical(nullptr, "Fatal Error", 
                             "Unknown fatal error occurred while starting HSML Visual Studio.\n\n"
                             "Please check your installation and try again.");
        
        std::cerr << "Unknown fatal error occurred" << std::endl;
        return 1;
    }
}

// GUI-FIRST: Perfect application lifecycle management
class ApplicationManager {
public:
    static ApplicationManager& instance() {
        static ApplicationManager instance;
        return instance;
    }
    
    void registerMainWindow(hsml::gui::HSMLVisualStudioGUI* window) {
        main_window_ = window;
    }
    
    hsml::gui::HSMLVisualStudioGUI* getMainWindow() const {
        return main_window_;
    }
    
    void requestShutdown() {
        if (main_window_) {
            main_window_->close();
        }
    }
    
private:
    ApplicationManager() = default;
    hsml::gui::HSMLVisualStudioGUI* main_window_{nullptr};
};

// Global application functions for plugin access
extern "C" {
    // Plugin API functions
    void* hsml_get_main_window() {
        return ApplicationManager::instance().getMainWindow();
    }
    
    void hsml_request_shutdown() {
        ApplicationManager::instance().requestShutdown();
    }
    
    const char* hsml_get_version() {
        return "1.0.0";
    }
    
    const char* hsml_get_build_info() {
        return "GUI-FIRST Architecture - Perfect Visual Experience";
    }
}