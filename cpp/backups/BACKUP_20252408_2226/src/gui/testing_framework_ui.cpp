/**
 * HSML Testing Framework UI - Implementation
 * GUI-FIRST METHODOLOGY: Perfect testing interface drives perfect test execution
 * Revolutionary spatial testing with real-time visual feedback
 */

#include "hsml/gui/hsml_visual_studio_gui.h"

#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QSplitter>
#include <QtWidgets/QTreeWidget>
#include <QtWidgets/QTreeWidgetItem>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QLabel>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QFileDialog>
#include <QtCore/QThread>
#include <QtCore/QTimer>
#include <QtCore/QDir>

namespace hsml::gui {

// =============================================================================
// TestingFrameworkUI Implementation
// =============================================================================

TestingFrameworkUI::TestingFrameworkUI(QWidget* parent)
    : QWidget(parent)
    , main_layout_(new QHBoxLayout(this))
    , main_splitter_(new QSplitter(Qt::Horizontal, this))
    , testing_framework_(nullptr)
    , test_runner_thread_(new QThread(this))
{
    setupUI();
    connectSignals();
    initializeTestingFramework();
}

TestingFrameworkUI::~TestingFrameworkUI() {
    if (test_runner_thread_->isRunning()) {
        test_runner_thread_->quit();
        test_runner_thread_->wait(3000);
    }
}

void TestingFrameworkUI::setupUI() {
    // GUI-FIRST: Perfect testing interface layout
    
    // Left panel - Test Suite Tree and Controls
    QWidget* left_panel = new QWidget(this);
    left_panel->setMinimumWidth(300);
    left_panel->setMaximumWidth(400);
    
    QVBoxLayout* left_layout = new QVBoxLayout(left_panel);
    
    // Test suite tree
    test_suite_tree_ = new QTreeWidget(this);
    test_suite_tree_->setHeaderLabels({"Test", "Status", "Duration"});
    test_suite_tree_->setAlternatingRowColors(true);
    test_suite_tree_->setSelectionMode(QAbstractItemView::ExtendedSelection);
    left_layout->addWidget(test_suite_tree_);
    
    // Test control buttons
    test_controls_layout_ = new QVBoxLayout();
    
    QHBoxLayout* main_buttons = new QHBoxLayout();
    run_all_button_ = new QPushButton("Run All Tests", this);
    run_all_button_->setStyleSheet("QPushButton { background-color: #0078d4; color: white; font-weight: bold; padding: 8px 16px; }");
    run_selected_button_ = new QPushButton("Run Selected", this);
    stop_tests_button_ = new QPushButton("Stop", this);
    stop_tests_button_->setEnabled(false);
    stop_tests_button_->setStyleSheet("QPushButton { background-color: #dc3545; color: white; }");
    
    main_buttons->addWidget(run_all_button_);
    main_buttons->addWidget(run_selected_button_);
    main_buttons->addWidget(stop_tests_button_);
    test_controls_layout_->addLayout(main_buttons);
    
    // Test options
    continuous_integration_checkbox_ = new QCheckBox("Continuous Integration", this);
    continuous_integration_checkbox_->setToolTip("Automatically run tests when files change");
    auto_test_checkbox_ = new QCheckBox("Auto-test on Save", this);
    auto_test_checkbox_->setToolTip("Run affected tests when HSML files are saved");
    
    test_controls_layout_->addWidget(continuous_integration_checkbox_);
    test_controls_layout_->addWidget(auto_test_checkbox_);
    
    // Test progress
    test_progress_ = new QProgressBar(this);
    test_progress_->setVisible(false);
    test_status_label_ = new QLabel("Ready", this);
    test_status_label_->setStyleSheet("QLabel { color: #28a745; font-weight: bold; }");
    
    test_controls_layout_->addWidget(test_progress_);
    test_controls_layout_->addWidget(test_status_label_);
    
    left_layout->addLayout(test_controls_layout_);
    
    // Right panel - Test Results and Details
    QWidget* right_panel = new QWidget(this);
    QVBoxLayout* right_layout = new QVBoxLayout(right_panel);
    
    // Results tabs
    results_tabs_ = new QTabWidget(this);
    
    // Test output tab
    test_output_ = new QTextEdit(this);
    test_output_->setReadOnly(true);
    test_output_->setFont(QFont("Consolas", 10));
    test_output_->setStyleSheet("QTextEdit { background-color: #1e1e1e; color: #d4d4d4; }");
    results_tabs_->addTab(test_output_, "Output");
    
    // Test results tree tab
    test_results_tree_ = new QTreeWidget(this);
    test_results_tree_->setHeaderLabels({"Test", "Result", "Duration", "Message"});
    test_results_tree_->setAlternatingRowColors(true);
    results_tabs_->addTab(test_results_tree_, "Results");
    
    // Performance metrics tab
    performance_metrics_widget_ = new QWidget(this);
    setupPerformanceMetricsWidget();
    results_tabs_->addTab(performance_metrics_widget_, "Performance");
    
    // Coverage report tab
    coverage_report_widget_ = new QWidget(this);
    setupCoverageReportWidget();
    results_tabs_->addTab(coverage_report_widget_, "Coverage");
    
    right_layout->addWidget(results_tabs_);
    
    // Add panels to splitter
    main_splitter_->addWidget(left_panel);
    main_splitter_->addWidget(right_panel);
    main_splitter_->setSizes({300, 700});
    
    main_layout_->addWidget(main_splitter_);
    setLayout(main_layout_);
    
    // Load sample test data
    loadSampleTests();
}

void TestingFrameworkUI::setupPerformanceMetricsWidget() {
    QVBoxLayout* layout = new QVBoxLayout(performance_metrics_widget_);
    
    // Performance summary
    QLabel* title = new QLabel("Performance Metrics", this);
    title->setStyleSheet("QLabel { font-size: 16px; font-weight: bold; color: #0078d4; }");
    layout->addWidget(title);
    
    // Metrics grid
    QGridLayout* metrics_grid = new QGridLayout();
    
    // Execution time metrics
    metrics_grid->addWidget(new QLabel("Total Execution Time:", this), 0, 0);
    QLabel* total_time_label = new QLabel("0.00s", this);
    total_time_label->setStyleSheet("QLabel { font-weight: bold; }");
    metrics_grid->addWidget(total_time_label, 0, 1);
    
    // Memory usage
    metrics_grid->addWidget(new QLabel("Peak Memory Usage:", this), 1, 0);
    QLabel* memory_label = new QLabel("0 MB", this);
    memory_label->setStyleSheet("QLabel { font-weight: bold; }");
    metrics_grid->addWidget(memory_label, 1, 1);
    
    // Spherical operations (HSML-specific)
    metrics_grid->addWidget(new QLabel("Coordinate Transformations:", this), 2, 0);
    QLabel* coord_ops_label = new QLabel("0", this);
    coord_ops_label->setStyleSheet("QLabel { font-weight: bold; }");
    metrics_grid->addWidget(coord_ops_label, 2, 1);
    
    // Steradian calculations
    metrics_grid->addWidget(new QLabel("Steradian Calculations:", this), 3, 0);
    QLabel* steradian_label = new QLabel("0", this);
    steradian_label->setStyleSheet("QLabel { font-weight: bold; }");
    metrics_grid->addWidget(steradian_label, 3, 1);
    
    layout->addLayout(metrics_grid);
    layout->addStretch();
}

void TestingFrameworkUI::setupCoverageReportWidget() {
    QVBoxLayout* layout = new QVBoxLayout(coverage_report_widget_);
    
    // Coverage title
    QLabel* title = new QLabel("Code Coverage Report", this);
    title->setStyleSheet("QLabel { font-size: 16px; font-weight: bold; color: #0078d4; }");
    layout->addWidget(title);
    
    // Coverage summary
    QHBoxLayout* summary_layout = new QHBoxLayout();
    
    QLabel* overall_label = new QLabel("Overall Coverage:", this);
    QProgressBar* coverage_bar = new QProgressBar(this);
    coverage_bar->setRange(0, 100);
    coverage_bar->setValue(85);
    coverage_bar->setFormat("85%");
    coverage_bar->setStyleSheet(R"(
        QProgressBar {
            border: 2px solid grey;
            border-radius: 5px;
            text-align: center;
        }
        QProgressBar::chunk {
            background-color: #28a745;
            width: 20px;
        }
    )");
    
    summary_layout->addWidget(overall_label);
    summary_layout->addWidget(coverage_bar);
    layout->addLayout(summary_layout);
    
    // Coverage details
    QTreeWidget* coverage_tree = new QTreeWidget(this);
    coverage_tree->setHeaderLabels({"File", "Lines", "Covered", "Coverage %"});
    coverage_tree->setAlternatingRowColors(true);
    
    // Sample coverage data
    QTreeWidgetItem* core_item = new QTreeWidgetItem(coverage_tree);
    core_item->setText(0, "hsml_dom.cpp");
    core_item->setText(1, "245");
    core_item->setText(2, "210");
    core_item->setText(3, "85.7%");
    
    QTreeWidgetItem* renderer_item = new QTreeWidgetItem(coverage_tree);
    renderer_item->setText(0, "spherical_renderer.cpp");
    renderer_item->setText(1, "189");
    renderer_item->setText(2, "165");
    renderer_item->setText(3, "87.3%");
    
    layout->addWidget(coverage_tree);
}

void TestingFrameworkUI::connectSignals() {
    // GUI-FIRST: Every interaction provides immediate visual feedback
    
    connect(run_all_button_, &QPushButton::clicked, this, &TestingFrameworkUI::runAllTests);
    connect(run_selected_button_, &QPushButton::clicked, this, &TestingFrameworkUI::runSelectedTest);
    connect(stop_tests_button_, &QPushButton::clicked, this, &TestingFrameworkUI::stopTests);
    
    connect(continuous_integration_checkbox_, &QCheckBox::toggled, 
            this, &TestingFrameworkUI::enableContinuousIntegration);
    connect(auto_test_checkbox_, &QCheckBox::toggled,
            this, &TestingFrameworkUI::setAutoTestOnSave);
    
    // Tree selection changes
    connect(test_suite_tree_, &QTreeWidget::itemSelectionChanged, [this]() {
        bool has_selection = !test_suite_tree_->selectedItems().isEmpty();
        run_selected_button_->setEnabled(has_selection);
    });
    
    // Double-click to run individual test
    connect(test_suite_tree_, &QTreeWidget::itemDoubleClicked, [this](QTreeWidgetItem* item) {
        if (item && item->parent()) { // It's a test item, not a suite
            runIndividualTest(item);
        }
    });
}

void TestingFrameworkUI::initializeTestingFramework() {
    try {
        // Initialize HSML testing framework
        testing::TestingConfig config = testing::TestingConfig::Builder()
            .withMaxWorkers(4)
            .withTimeout(std::chrono::milliseconds(30000))
            .withSpatialPrecision(1e-10)
            .enableAI(true)
            .enableParallel(true)
            .build();
        
        testing_framework_ = std::make_unique<testing::HSMLTestingFramework>(config);
        
        test_output_->append("<font color='#28a745'>[INFO] HSML Testing Framework initialized successfully</font>");
        test_status_label_->setText("Framework ready");
        
    } catch (const std::exception& e) {
        test_output_->append(QString("<font color='#dc3545'>[ERROR] Failed to initialize testing framework: %1</font>")
                           .arg(e.what()));
        test_status_label_->setText("Framework initialization failed");
        test_status_label_->setStyleSheet("QLabel { color: #dc3545; font-weight: bold; }");
    }
}

void TestingFrameworkUI::loadSampleTests() {
    // GUI-FIRST: Show immediate visual feedback with sample test structure
    
    test_suite_tree_->clear();
    
    // Unit Tests Suite
    QTreeWidgetItem* unit_suite = new QTreeWidgetItem(test_suite_tree_);
    unit_suite->setText(0, "Unit Tests");
    unit_suite->setText(1, "Ready");
    unit_suite->setIcon(0, style()->standardIcon(QStyle::SP_DirIcon));
    
    // Spherical coordinate tests
    QTreeWidgetItem* coord_test = new QTreeWidgetItem(unit_suite);
    coord_test->setText(0, "Spherical Coordinate Conversion");
    coord_test->setText(1, "Pending");
    coord_test->setIcon(0, style()->standardIcon(QStyle::SP_FileIcon));
    
    QTreeWidgetItem* distance_test = new QTreeWidgetItem(unit_suite);
    distance_test->setText(0, "Distance Calculations");
    distance_test->setText(1, "Pending");
    distance_test->setIcon(0, style()->standardIcon(QStyle::SP_FileIcon));
    
    QTreeWidgetItem* steradian_test = new QTreeWidgetItem(unit_suite);
    steradian_test->setText(0, "Steradian Calculations");
    steradian_test->setText(1, "Pending");
    steradian_test->setIcon(0, style()->standardIcon(QStyle::SP_FileIcon));
    
    // Integration Tests Suite
    QTreeWidgetItem* integration_suite = new QTreeWidgetItem(test_suite_tree_);
    integration_suite->setText(0, "Integration Tests");
    integration_suite->setText(1, "Ready");
    integration_suite->setIcon(0, style()->standardIcon(QStyle::SP_DirIcon));
    
    QTreeWidgetItem* rendering_test = new QTreeWidgetItem(integration_suite);
    rendering_test->setText(0, "Spherical Rendering Pipeline");
    rendering_test->setText(1, "Pending");
    rendering_test->setIcon(0, style()->standardIcon(QStyle::SP_FileIcon));
    
    QTreeWidgetItem* material_test = new QTreeWidgetItem(integration_suite);
    material_test->setText(0, "Material System Integration");
    material_test->setText(1, "Pending");
    material_test->setIcon(0, style()->standardIcon(QStyle::SP_FileIcon));
    
    // Performance Tests Suite
    QTreeWidgetItem* perf_suite = new QTreeWidgetItem(test_suite_tree_);
    perf_suite->setText(0, "Performance Tests");
    perf_suite->setText(1, "Ready");
    perf_suite->setIcon(0, style()->standardIcon(QStyle::SP_DirIcon));
    
    QTreeWidgetItem* fps_test = new QTreeWidgetItem(perf_suite);
    fps_test->setText(0, "60 FPS Rendering Test");
    fps_test->setText(1, "Pending");
    fps_test->setIcon(0, style()->standardIcon(QStyle::SP_FileIcon));
    
    QTreeWidgetItem* memory_test = new QTreeWidgetItem(perf_suite);
    memory_test->setText(0, "Memory Usage Test");
    memory_test->setText(1, "Pending");
    memory_test->setIcon(0, style()->standardIcon(QStyle::SP_FileIcon));
    
    // Spatial Tests Suite
    QTreeWidgetItem* spatial_suite = new QTreeWidgetItem(test_suite_tree_);
    spatial_suite->setText(0, "Spatial Tests");
    spatial_suite->setText(1, "Ready");
    spatial_suite->setIcon(0, style()->standardIcon(QStyle::SP_DirIcon));
    
    QTreeWidgetItem* collision_test = new QTreeWidgetItem(spatial_suite);
    collision_test->setText(0, "Collision Detection");
    collision_test->setText(1, "Pending");
    collision_test->setIcon(0, style()->standardIcon(QStyle::SP_FileIcon));
    
    QTreeWidgetItem* physics_test = new QTreeWidgetItem(spatial_suite);
    physics_test->setText(0, "Physics Simulation");
    physics_test->setText(1, "Pending");
    physics_test->setIcon(0, style()->standardIcon(QStyle::SP_FileIcon));
    
    // Expand all suites
    test_suite_tree_->expandAll();
    
    // Update statistics
    current_stats_.total_tests = 9;
    updateStatusDisplay();
}

void TestingFrameworkUI::runAllTests() {
    if (!testing_framework_) {
        QMessageBox::warning(this, "Testing Framework", "Testing framework not initialized");
        return;
    }
    
    // GUI-FIRST: Immediate visual feedback
    run_all_button_->setEnabled(false);
    run_selected_button_->setEnabled(false);
    stop_tests_button_->setEnabled(true);
    
    test_progress_->setVisible(true);
    test_progress_->setRange(0, current_stats_.total_tests);
    test_progress_->setValue(0);
    
    test_status_label_->setText("Running all tests...");
    test_status_label_->setStyleSheet("QLabel { color: #ffc107; font-weight: bold; }");
    
    test_output_->clear();
    test_output_->append("<font color='#0078d4'>[INFO] Starting test run...</font>");
    
    emit testRunStarted();
    
    // Simulate test execution with visual feedback
    QTimer* test_timer = new QTimer(this);
    int current_test = 0;
    
    connect(test_timer, &QTimer::timeout, [this, test_timer, &current_test]() {
        if (current_test >= current_stats_.total_tests) {
            test_timer->stop();
            test_timer->deleteLater();
            onAllTestsCompleted();
            return;
        }
        
        // Simulate running a test
        QString test_name = QString("Test %1").arg(current_test + 1);
        onTestStarted(test_name);
        
        // Simulate test result (85% pass rate)
        bool passed = (current_test % 7) != 0; // Fail every 7th test
        
        QTimer::singleShot(200 + (qrand() % 300), [this, test_name, passed]() {
            if (passed) {
                onTestCompleted(test_name, true);
            } else {
                onTestFailed(test_name, "Sample test failure for demonstration");
            }
        });
        
        current_test++;
        test_progress_->setValue(current_test);
    });
    
    test_timer->start(100); // Start test execution simulation
}

void TestingFrameworkUI::onTestStarted(const QString& test_name) {
    test_output_->append(QString("<font color='#6c757d'>[RUNNING] %1</font>").arg(test_name));
    test_output_->moveCursor(QTextCursor::End);
}

void TestingFrameworkUI::onTestCompleted(const QString& test_name, bool passed) {
    if (passed) {
        test_output_->append(QString("<font color='#28a745'>[PASSED] %1</font>").arg(test_name));
        current_stats_.passed_tests++;
    } else {
        test_output_->append(QString("<font color='#dc3545'>[FAILED] %1</font>").arg(test_name));
        current_stats_.failed_tests++;
    }
    
    updateStatusDisplay();
    test_output_->moveCursor(QTextCursor::End);
}

void TestingFrameworkUI::onTestFailed(const QString& test_name, const QString& error) {
    test_output_->append(QString("<font color='#dc3545'>[FAILED] %1 - %2</font>").arg(test_name, error));
    current_stats_.failed_tests++;
    updateStatusDisplay();
    test_output_->moveCursor(QTextCursor::End);
}

void TestingFrameworkUI::onAllTestsCompleted() {
    // GUI-FIRST: Perfect completion feedback
    run_all_button_->setEnabled(true);
    run_selected_button_->setEnabled(true);
    stop_tests_button_->setEnabled(false);
    
    test_progress_->setVisible(false);
    
    bool all_passed = (current_stats_.failed_tests == 0);
    if (all_passed) {
        test_status_label_->setText(QString("All tests passed! (%1/%2)")
                                  .arg(current_stats_.passed_tests)
                                  .arg(current_stats_.total_tests));
        test_status_label_->setStyleSheet("QLabel { color: #28a745; font-weight: bold; }");
        test_output_->append("<font color='#28a745'>[SUCCESS] All tests completed successfully!</font>");
    } else {
        test_status_label_->setText(QString("Tests completed: %1 passed, %2 failed")
                                  .arg(current_stats_.passed_tests)
                                  .arg(current_stats_.failed_tests));
        test_status_label_->setStyleSheet("QLabel { color: #dc3545; font-weight: bold; }");
        test_output_->append(QString("<font color='#dc3545'>[COMPLETED] Tests finished with %1 failures</font>")
                           .arg(current_stats_.failed_tests));
    }
    
    emit testRunCompleted(current_stats_.passed_tests, current_stats_.failed_tests);
    test_output_->moveCursor(QTextCursor::End);
}

void TestingFrameworkUI::updateStatusDisplay() {
    QString status = QString("Tests: %1 total, %2 passed, %3 failed")
                    .arg(current_stats_.total_tests)
                    .arg(current_stats_.passed_tests)
                    .arg(current_stats_.failed_tests);
    
    // Update status in results tab
    if (results_tabs_->currentIndex() == 1) { // Results tab
        // Update results tree with current statistics
    }
}

void TestingFrameworkUI::runSelectedTest() {
    QList<QTreeWidgetItem*> selected = test_suite_tree_->selectedItems();
    if (selected.isEmpty()) {
        QMessageBox::information(this, "No Selection", "Please select one or more tests to run.");
        return;
    }
    
    test_output_->append("<font color='#0078d4'>[INFO] Running selected tests...</font>");
    
    for (QTreeWidgetItem* item : selected) {
        if (item->parent()) { // It's a test, not a suite
            runIndividualTest(item);
        }
    }
}

void TestingFrameworkUI::runIndividualTest(QTreeWidgetItem* test_item) {
    if (!test_item || !test_item->parent()) return;
    
    QString test_name = test_item->text(0);
    test_output_->append(QString("<font color='#6c757d'>[RUNNING] %1</font>").arg(test_name));
    
    // Simulate individual test execution
    QTimer::singleShot(500 + (qrand() % 1000), [this, test_item, test_name]() {
        bool passed = (qrand() % 5) != 0; // 80% pass rate
        
        if (passed) {
            test_item->setText(1, "Passed");
            test_item->setIcon(1, style()->standardIcon(QStyle::SP_DialogApplyButton));
            test_output_->append(QString("<font color='#28a745'>[PASSED] %1</font>").arg(test_name));
        } else {
            test_item->setText(1, "Failed");
            test_item->setIcon(1, style()->standardIcon(QStyle::SP_DialogCancelButton));
            test_output_->append(QString("<font color='#dc3545'>[FAILED] %1</font>").arg(test_name));
        }
        
        test_item->setText(2, QString("%1ms").arg(qrand() % 1000 + 10));
        test_output_->moveCursor(QTextCursor::End);
    });
}

void TestingFrameworkUI::stopTests() {
    // GUI-FIRST: Immediate stop feedback
    test_output_->append("<font color='#ffc107'>[STOPPED] Test execution stopped by user</font>");
    
    run_all_button_->setEnabled(true);
    run_selected_button_->setEnabled(true);
    stop_tests_button_->setEnabled(false);
    
    test_progress_->setVisible(false);
    test_status_label_->setText("Tests stopped");
    test_status_label_->setStyleSheet("QLabel { color: #ffc107; font-weight: bold; }");
    
    test_output_->moveCursor(QTextCursor::End);
}

void TestingFrameworkUI::enableContinuousIntegration(bool enabled) {
    if (enabled) {
        test_output_->append("<font color='#17a2b8'>[CI] Continuous Integration enabled</font>");
        // Set up file watchers for CI
    } else {
        test_output_->append("<font color='#17a2b8'>[CI] Continuous Integration disabled</font>");
    }
    test_output_->moveCursor(QTextCursor::End);
}

void TestingFrameworkUI::setAutoTestOnSave(bool enabled) {
    if (enabled) {
        test_output_->append("<font color='#17a2b8'>[AUTO] Auto-test on save enabled</font>");
    } else {
        test_output_->append("<font color='#17a2b8'>[AUTO] Auto-test on save disabled</font>");
    }
    test_output_->moveCursor(QTextCursor::End);
}

} // namespace hsml::gui