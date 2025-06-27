#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

struct Candle {
    double close;
};

struct TradeResult {
    double success_rate;
    double avg_return;
    int total_trades;
    vector<int> signals; 
};

// rsi
double calculate_rsi(const vector<double>& closes, int index, int period = 14) {
    if (index < period) return 50.0;
    double gain = 0.0, loss = 0.0;
    for (int i = index - period + 1; i <= index; ++i) {
        double diff = closes[i] - closes[i - 1];
        if (diff > 0) gain += diff;
        else loss -= diff;
    }
    double rs = gain / (loss == 0 ? 1 : loss);
    return 100.0 - (100.0 / (1 + rs));
}

TradeResult rsi_strategy(const vector<Candle>& candles, double threshold) {
    vector<double> closes;
    for (auto c : candles) closes.push_back(c.close);

    int trades = 0, profit_trades = 0;
    double total_return = 0.0;
    bool in_position = false;
    double entry_price = 0;
    vector<int> signals(closes.size(), 0);

    for (int i = 15; i < closes.size(); ++i) {
        double rsi = calculate_rsi(closes, i);
        if (!in_position && rsi < 30) {
            entry_price = closes[i];
            in_position = true;
            signals[i] = 1;  // buy
        } else if (in_position && rsi > 70) {
            double ret = (closes[i] - entry_price) / entry_price;
            total_return += ret;
            if (ret > threshold) profit_trades++;
            trades++;
            in_position = false;
            signals[i] = -1; // sell
        }
    }
    double avg = trades ? total_return / trades * 100 : 0;
    double rate = trades ? profit_trades * 100.0 / trades : 0;
    return {rate, avg, trades, signals};
}

// MACD
double ema(const vector<double>& prices, int current, int period) {
    double k = 2.0 / (period + 1);
    double ema = prices[current - period];
    for (int i = current - period + 1; i <= current; ++i)
        ema = prices[i] * k + ema * (1 - k);
    return ema;
}

TradeResult macd_strategy(const vector<Candle>& candles, double threshold) {
    vector<double> closes;
    for (auto c : candles) closes.push_back(c.close);

    vector<double> macd;
    for (int i = 26; i < closes.size(); ++i) {
        double ema12 = ema(closes, i, 12);
        double ema26 = ema(closes, i, 26);
        macd.push_back(ema12 - ema26);
    }

    vector<double> signal;
    for (int i = 9; i < macd.size(); ++i) {
        double sum = 0;
        for (int j = i - 9 + 1; j <= i; ++j) sum += macd[j];
        signal.push_back(sum / 9);
    }

    int trades = 0, profit_trades = 0;
    double total_return = 0.0;
    bool in_position = false;
    double entry_price = 0;
    vector<int> signals(closes.size(), 0);

    for (int i = 1; i < signal.size(); ++i) {
        double m = macd[i+8];
        double s = signal[i];
        double pm = macd[i+7];
        double ps = signal[i-1];
        int price_index = i + 26 + 8;

        if (price_index >= closes.size()) break;

        if (!in_position && pm < ps && m > s) {
            entry_price = closes[price_index];
            in_position = true;
            signals[price_index] = 1;
        } else if (in_position && pm > ps && m < s) {
            double ret = (closes[price_index] - entry_price) / entry_price;
            total_return += ret;
            if (ret > threshold) profit_trades++;
            trades++;
            in_position = false;
            signals[price_index] = -1; 
        }
    }
    double avg = trades ? total_return / trades * 100 : 0;
    double rate = trades ? profit_trades * 100.0 / trades : 0;
    return {rate, avg, trades, signals};
}

// SMA
double sma(const vector<double>& closes, int end, int period) {
    if (end < period) return 0;
    double sum = 0;
    for (int i = end - period + 1; i <= end; ++i) sum += closes[i];
    return sum / period;
}

TradeResult sma_crossover_strategy(const vector<Candle>& candles, double threshold) {
    vector<double> closes;
    for (auto c : candles) closes.push_back(c.close);

    int trades = 0, profit_trades = 0;
    double total_return = 0;
    bool in_position = false;
    double entry_price = 0;
    vector<int> signals(closes.size(), 0);

    for (int i = 50; i < closes.size(); ++i) {
        double sma20 = sma(closes, i, 20);
        double sma50 = sma(closes, i, 50);
        double prev_sma20 = sma(closes, i - 1, 20);
        double prev_sma50 = sma(closes, i - 1, 50);

        if (!in_position && prev_sma20 < prev_sma50 && sma20 > sma50) {
            entry_price = closes[i];
            in_position = true;
            signals[i] = 1;
        } else if (in_position && prev_sma20 > prev_sma50 && sma20 < sma50) {
            double ret = (closes[i] - entry_price) / entry_price;
            total_return += ret;
            if (ret > threshold) profit_trades++;
            trades++;
            in_position = false;
            signals[i] = -1;
        }
    }
    double avg = trades ? total_return / trades * 100 : 0;
    double rate = trades ? profit_trades * 100.0 / trades : 0;
    return {rate, avg, trades, signals};
}
